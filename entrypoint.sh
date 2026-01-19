#!/usr/bin/env bash
set -euo pipefail

# Ensure writable dirs exist (useful also when run by Nextflow container engine)
mkdir -p "${NXF_HOME:-/work/.nextflow}" "${NXF_WORK:-/work/work}" "${NXF_TEMP:-/tmp}" || true

# Ensure Java uses a writable temp dir
export NXF_OPTS="${NXF_OPTS:-} -Djava.io.tmpdir=${NXF_TEMP:-/tmp}"

# Pass-through: execute exactly what is provided as CMD/args
exec "$@"