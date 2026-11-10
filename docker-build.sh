#!/usr/bin/env bash
# Build the GLIMPSE2 Docker image from a commit or tag of this repository.
#
# Usage: ./docker-build.sh [ref] [-t tag] [--push repository]
#
#   ref                 commit, tag or branch to build (default HEAD); needs this Dockerfile (2.1+)
#   -t tag              image tag (default glimpse:<version> for tag v<version>, else
#                       glimpse:<version>-<commit>)
#   --push repository   build linux/amd64 and linux/arm64 with buildx and push to repository
#
# The ref is checked out in a temporary worktree, so uncommitted changes are not built.
set -euo pipefail

usage() { sed -n '4,10s/^# \{0,1\}//p' "$0"; exit "${1:-0}"; }

ref=HEAD tag="" push=""
while [ $# -gt 0 ]; do
  case $1 in
    -t) tag=${2:?-t needs a tag}; shift 2 ;;
    --push) push=${2:?--push needs a repository}; shift 2 ;;
    -h|--help) usage ;;
    -*) echo "unknown option $1" >&2; usage 1 ;;
    *) ref=$1; shift ;;
  esac
done

repo=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
src=$(mktemp -d)
# git worktree remove refuses worktrees with submodules, and submodule deinit would also
# deinitialise SIMDe in the main checkout (shared config), so delete and prune instead.
cleanup() { rm -rf "$src"; git -C "$repo" worktree prune; }
trap cleanup EXIT

git -C "$repo" worktree add --quiet --detach "$src" "$ref"
simde=$(git -C "$repo" rev-parse --path-format=absolute --git-common-dir)/modules/third_party/simde
if [ -d "$simde" ]; then
  git -C "$src" submodule update --quiet --init --reference "$simde"
else
  git -C "$src" submodule update --quiet --init
fi
grep -q '^ARG GLIMPSE_COMMIT' "$src/Dockerfile" || { echo "$ref predates this Dockerfile" >&2; exit 1; }

commit=$(git -C "$src" rev-parse --short HEAD)
date=$(git -C "$src" log -1 --format=%cd --date=short)
version=$(sed -n 's/^#define PHASE_VERSION "\(.*\)"$/\1/p' "$src/versions/versions.h")
if [ -z "$tag" ]; then
  name=${push:-glimpse}
  if git -C "$src" tag --points-at HEAD | grep -qx "v$version"; then tag=$name:$version; else tag=$name:$version-$commit; fi
fi

echo "Building $tag from $ref ($commit, $date)"
args=(-t "$tag" --build-arg GLIMPSE_COMMIT="$commit" --build-arg GLIMPSE_DATE="$date")
if [ -n "$push" ]; then
  docker buildx build --platform linux/amd64,linux/arm64 --push "${args[@]}" "$src"
else
  docker build "${args[@]}" "$src"
fi
