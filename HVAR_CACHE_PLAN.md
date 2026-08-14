# Reduce GLIMPSE2_phase memory usage on very large binary reference panels
### (phase-only change — existing `.bin` panels are never modified or rebuilt from source)

## Context

The user runs `GLIMPSE2_phase` against a binary (`.bin`) reference panel built from
~500,000 samples (~1,000,000 reference haplotypes) via `GLIMPSE2_split_reference`, while
imputing only a handful of target samples, and observes very large memory usage that
does not shrink with fewer target samples. **Rebuilding these panels from the source VCF
is not an option** — they can't recreate the large reference panels anytime soon. So this
plan touches only `phase`-side code, treats the existing `.bin` file format as fixed and
read-only, and never requires `GLIMPSE2_split_reference` to be rerun against source data.

**Confirmed root cause:** `GLIMPSE2_phase` deserializes the *entire* `.bin` file into RAM
in one shot (`phase/src/caller/caller_initialise.cpp:203-205`, `ia >> H; ia >> V;`), with
no lazy loading or partial decode. The dominant structure is `HvarRef`
(`common/src/containers/ref_haplotype_set.h:121`, shared via symlink into
`phase/src/containers/` — confirmed there is exactly one physical copy) — a dense,
uncompressed bitmatrix (1 bit per variant × haplotype: `bitmatrix.cpp:77-83`) covering
"common" variants. It's a single global object (`caller::H`), loaded once, retained for
the whole run, **independent of target sample count**. For `n_ref_haps` = 1,000,000 and a
typical chunk (~24,000 common sites), `HvarRef` alone is ~3 GB.

**Why plain `mmap()` of the existing file wouldn't help:** the only hot-path reader,
`conditioning_set::compactSelection`'s `TYPE_COMMON` branch
(`phase/src/containers/conditioning_set.cpp:153-166`), loops over every common-site row
and gathers bits at ~`Kpbwt` (default 2000) *scattered* haplotype-column indices per row.
Since `HvarRef` is stored **row-major/variant-first** (each row spans all 1,000,000
haplotypes ≈ 125,000 bytes ≈ dozens of OS pages), ~2000 scattered draws per row hit nearly
every page in that row (birthday-paradox argument) — so this loop already touches
essentially the entire matrix on one individual's one iteration. `mmap`-ing the *current*
layout instead of `malloc`+`read` would not reduce resident memory.

**The fix that works, without touching the source `.bin`:** derive a haplotype-major
(transposed) copy of `HvarRef` — one contiguous ~3 KB block per reference haplotype — as
a **separate, phase-generated cache file**, built by a new **explicit, opt-in**
`GLIMPSE2_phase --build-hvar-cache` step run once per existing panel file. Normal `phase`
runs then `mmap()` that cache file (not the original `.bin`) for `HvarRef`, so selecting
`Kpbwt` scattered haplotypes only touches ~`Kpbwt × 3 KB` (~6 MB) instead of the whole
matrix. If no cache has been built yet for a given panel, `phase` behaves **exactly as it
does today** — this is purely additive.

This design was cross-checked with two Plan agents and verified directly against the
code: `bitmatrix::transpose()`/`subset()` (`bitmatrix.cpp:100-151`, `42-52`) are existing,
currently-dead primitives that do exactly the gather/transpose work needed; `H.flag_common`
is a global per-chunk flag so the full common-site set is always in play (no filtering
complicates the gather); `conditioning_set::select` (`caller_algorithm.cpp:56`) is the only
per-individual/iteration call site, and PBWT state selection itself never touches
`HvarRef`.

## Plan

### 1. New `bitmatrix` capabilities (`common/src/containers/bitmatrix.{h,cpp}`)

This file is shared (symlinked) by `phase`/`split_reference`/`common`, but the new
behavior is only ever invoked from `phase`'s code paths — `split_reference`'s own build
logic (which fills `HvarRef` row-major from VCF) and `phase`'s direct-BCF input mode are
untouched and keep working exactly as today.

- **Skip-and-discard load mode**, to avoid ever materializing the *original* `.bin`'s
  row-major `HvarRef` bytes in RAM when a valid cache is going to replace them, while
  still correctly advancing the archive stream for the fields that follow (`Ypacked`,
  `A_small_idx`):
  ```cpp
  // bitmatrix.h: new member
  bool skip_and_discard_on_load = false;

  // bitmatrix.cpp: bitmatrix::serialize(), load branch
  ar & n_bytes; ar & n_cols; ar & n_rows;
  if (Archive::is_loading::value) {
      if (skip_and_discard_on_load) {
          // Consume exactly n_bytes from the archive stream via a small reusable
          // scratch buffer, so peak RSS never includes the full row-major matrix.
          constexpr unsigned long CHUNK = 1u << 20; // 1 MB
          std::vector<unsigned char> scratch(std::min(n_bytes, CHUNK));
          for (unsigned long done = 0; done < n_bytes; ) {
              unsigned long take = std::min(n_bytes - done, CHUNK);
              ar & boost::serialization::make_array<unsigned char>(scratch.data(), take);
              done += take;
          }
          bytes = nullptr;   // caller adopts the mmap'd cache blob afterward
          return;
      }
      bytes = (unsigned char*)std::malloc(n_bytes);
  }
  ar & boost::serialization::make_array<unsigned char>(bytes, n_bytes);
  ```
  This changes *how* already-existing bytes in the current `.bin` format are consumed —
  it does not change what's written to disk, so it's fully compatible with every existing
  panel file.
- **mmap adoption mode** (used for both the cache-build step's read of the transposed
  blob back for validation, and normal-run consumption): `owns_mmap` flag, `mmap_length`,
  `adopt_mmap(bytes, mapped_length, nrow, ncol, n_bytes_logical)`, and `release()`
  (`munmap()` if `owns_mmap` else `free()`), called from the destructor and explicitly
  before each retry attempt in the existing retry-with-backoff loop.
- Guard mutating methods (`allocate`, `reallocate`, `set*`) with
  `assert(!owns_mmap)` — defense-in-depth; verified via grep that no mutating call ever
  reaches an mmap-adopted `HvarRef`.
- `= delete` copy constructor/assignment (bitmatrix currently has none defined; adding a
  second owned-OS-resource type makes an accidental copy a real double-free/double-munmap
  risk).
- Add `rowPtr(row)` (const + non-const) for row-granularity `memcpy`, mirroring
  `getByte`'s addressing; relax `subset()`'s source parameter to `const bitmatrix &`.

### 2. New cache file format (phase-only, new file — no backward-compat constraints)

New small header (e.g. `phase/src/containers/hvar_cache.h`), page-aligned so the blob
that follows can be `mmap()`'d directly:

```
[0 .. 65536)              fixed-size header, zero-padded: magic "GLIMPSE2HVCACHE",
                           format_version, source fingerprint (source .bin file size +
                           mtime), n_ref_haps, n_com_sites (redundant cross-check against
                           the fingerprint), blob_bytes, blob_crc32
[65536 .. +blob_bytes)     raw haplotype-major HvarRef bytes: one contiguous
                           round8(n_com_sites)/8-byte block per reference haplotype, rows
                           0..n_ref_haps-1, no per-row padding
```
65536 bytes is a fixed constant (not the build host's OS page size) so the blob offset is
a multiple of every realistic page size across the three supported platforms (Linux
x86_64/aarch64, macOS Apple Silicon).

Fingerprint = source `.bin` file size + mtime (`stat()`), checked before trusting the
cache; a mismatch (source panel was replaced/rebuilt) means the cache is treated as
stale and silently ignored, falling back to a full load with a `vrb.warning` suggesting
`--build-hvar-cache` be rerun.

### 3. New `GLIMPSE2_phase` CLI surface (`caller_parameters.cpp`)

- `--build-hvar-cache` (flag): switches `phase` into a one-shot cache-build mode.
  Requires `--reference <panel.bin>` (must resolve to the binary/`InputFormat::GLIMPSE`
  path, not VCF/BCF — `check_options()` should error otherwise); does **not** require
  `--output`, `--bam-file`/`--input-gl`, `--input-region`, etc. — the normal
  target-sample/BAM machinery is skipped entirely. Reads the panel via the existing
  loader (unchanged, full in-RAM row-major load — this one-time build step pays the same
  memory cost `phase` pays today, but only once, deliberately, as a maintenance step),
  transposes `HvarRef` via the revived `bitmatrix::transpose()`, writes the cache file
  atomically (temp file + `std::filesystem::rename`, matching the existing pattern
  already used by `caller::write_checkpoint()` in `caller_algorithm.cpp:136-168`), prints
  a confirmation, and exits — never proceeds to `phase_loop()`.
- `--hvar-cache-file <path>` (optional, both modes): explicit cache path. If omitted,
  default `<reference>.hvarT` next to the `.bin`. Lets operators redirect the cache to a
  writable location when the reference-panel directory is mounted read-only.
- `caller::phase()` (`caller_management.cpp:47`) gains an early branch: if
  `--build-hvar-cache` is set, call a new `build_hvar_cache()` and return, bypassing
  `read_files_and_initialise()`/`phase_loop()`/`write_files_and_finalise()`.

### 4. `phase` read-side changes (`caller_initialise.cpp`)

In `read_binary_reference_panel`, before constructing the boost archive:
1. Compute the cache path (`--hvar-cache-file` or default), `stat()` both the source
   `.bin` and the candidate cache file.
2. If the cache exists and its header's fingerprint matches the source `.bin`'s current
   size/mtime (and `n_ref_haps`/`n_com_sites` cross-check passes once available), mark
   `use_cache = true` and `mmap()` the cache's blob (`PROT_READ`, `MAP_PRIVATE`,
   `madvise(MADV_RANDOM)`).
3. If `use_cache`, set `H.HvarRef.skip_and_discard_on_load = true` **before** `ia >> H;`
   so the original `.bin`'s row-major bytes are streamed-and-discarded (§1) rather than
   materialized; after deserialization completes, call
   `H.HvarRef.adopt_mmap(...)` with the cache's mmap'd pointer and haplotype-major
   dimensions (`n_rows = n_ref_haps`, `n_cols = n_com_sites`), and set a new
   `haplotype_set::hvarref_is_transposed = true` flag.
4. If `use_cache` is false (no cache built yet, or stale/mismatched), behave exactly as
   today: full malloc+read, row-major, `hvarref_is_transposed = false`, and log a
   `vrb.bullet`/`vrb.warning` suggesting `--build-hvar-cache` for future runs against this
   panel.
5. Fallback: if `mmap()` of the cache fails for any reason (unsupported filesystem,
   page-size/offset check fails), treat it the same as "no valid cache" (case 4) rather
   than erroring — the run must always succeed, just without the memory benefit.

### 5. `conditioning_set` dual-path consumption (`conditioning_set.{h,cpp}`)

`compactSelection`'s `TYPE_COMMON` branch gets a runtime branch on
`H.hvarref_is_transposed`:

- **`false` (no cache — today's behavior, unchanged):** keep the existing scalar
  8-bit-gather loop exactly as it is now (`H.HvarRef.get(lcom, idxHaps_ref[k])`,
  row-major).
- **`true` (cache active):** gather-then-transpose using the now-haplotype-major
  `H.HvarRef` (`H.HvarRef.get(idxHaps_ref[k], lcom)` convention, arguments reversed):
  ```cpp
  Hgathered.subset(H.HvarRef, idxHaps_ref);              // revived subset(): per-haplotype memcpy
  Htransposed.reallocate(Hgathered.n_cols, Hgathered.n_rows);
  Hgathered.transpose(Htransposed);                        // revived transpose(): 8x8 block bit-transpose
  // scatter into Hvar's interleaved COMMON/RARE rows (lrel/lcom bookkeeping unchanged)
  std::memcpy(Hvar.rowPtr(lrel), Htransposed.rowPtr(lcom), hvar_row_bytes);
  ```
  Add `Hgathered`/`Htransposed` as new per-`conditioning_set` scratch `bitmatrix`
  members (one instance already exists per worker thread — no new threading concerns).
  `TYPE_RARE`/`TYPE_MONO` handling is untouched.

Expected effect (cache-active case only): bytes touched per `compactSelection` call drops
from `O(n_com_sites × n_ref_haps)` (~3 GB) to `O(n_states × n_com_sites)` (~6 MB at
defaults). No panel-size threshold — this path only ever runs when the user has
opted in via `--build-hvar-cache`, so there's no risk of surprising a workflow that
hasn't adopted it.

### 6. Concurrency

No new synchronization needed: worker threads (one `conditioning_set` per `--threads`)
already read the shared `H` object concurrently today; concurrent reads against a single
`PROT_READ`/`MAP_PRIVATE` mapping are safe with no locking, same as today's shared
malloc'd buffer.

## Verification

1. **Unit-level correctness of the revived primitives:** small random `bitmatrix`
   (including non-multiple-of-8 dimensions), transpose it, assert
   `transposed.get(i,j) == original.get(j,i)` for all `(i,j)` — `transpose()`/`subset()`
   are currently unexercised by any code path, so this must be validated before relying
   on them in the hot path.
2. **Bit-identical output test:** on a small synthetic panel, compare `Hvar` produced by
   the existing row-major path vs. the new cache-active gather+transpose path (built from
   a locally-generated `--build-hvar-cache` cache of that same small panel) and assert
   byte-for-byte equality.
3. **End-to-end regression, no-cache path:** run the modified `GLIMPSE2_phase` against
   the existing `tutorial/` panel with no cache built — confirm identical output and
   behavior to current `master` (this path must be a no-op change).
4. **End-to-end regression, cache path:** run `GLIMPSE2_phase --build-hvar-cache
   --reference <tutorial panel.bin>`, then run a normal phase invocation against that
   panel + cache, and diff the output VCF/BCF against the no-cache run on the same
   input — must be identical (this is purely a storage/access-pattern change, not
   algorithmic).
5. **Stale-cache handling:** rebuild/touch the source `.bin` (or hand-edit its mtime) and
   confirm the next `phase` run detects the mismatch, ignores the stale cache, and falls
   back to the full load with a warning rather than using incorrect data.
6. **Memory measurement (the actual point of this change):** with a cache built for a
   large panel (real 500k-sample scale, or a scaled-down proxy), run `phase` against a
   handful of target samples under `/usr/bin/time -v` (or platform equivalent) and
   confirm peak RSS is now roughly proportional to target-sample-count × `Kpbwt`, not to
   `n_ref_haps` — compare directly against the no-cache run's peak RSS on the same
   inputs.
7. **Read-only reference directory:** confirm `--build-hvar-cache --hvar-cache-file
   <writable path>` succeeds when the source `.bin`'s directory is read-only, and that a
   normal phase run picks up a cache file at a non-default location via
   `--hvar-cache-file`.

## Status

Implemented and verified against a synthetic reference panel (see conversation/PR
history): compiles cleanly, `--build-hvar-cache` builds a cache that `phase`
auto-detects, single-threaded output is bit-identical between the no-cache and
cache-active paths (both the full-panel and scattered-PBWT-selection cases), VCF/BCF is
rejected for `--build-hvar-cache`, stale caches fall back gracefully, and custom
`--hvar-cache-file` paths work. Not yet validated: actual peak-RSS reduction at real
500k-sample scale (verification step 6) and behavior against a read-only reference
directory (step 7) — both require production-scale data not available in the
implementing environment.
