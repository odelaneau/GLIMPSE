zcat $1 | awk 'BEGIN { printf("["); } { printf("{\"x\":%.8f,\"y\":%.4f}", $3*100, $5) } END { printf("]"); }' | sed 's/}{/}, {/g'
