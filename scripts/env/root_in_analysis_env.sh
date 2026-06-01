#!/usr/bin/env bash

set -euo pipefail

ENVROOT="${HOME}/Desktop/analysis/env"
ROOUNFOLD_OLD="${HOME}/Desktop/ThesisAnalysis/external/RooUnfold"

if [[ ! -x "${ENVROOT}/bin/root" ]]; then
  echo "Missing ROOT binary at ${ENVROOT}/bin/root" >&2
  exit 1
fi

if command -v xcrun >/dev/null 2>&1; then
  if sdkroot="$(xcrun --sdk macosx --show-sdk-path 2>/dev/null)"; then
    export SDKROOT="${sdkroot}"
  fi
fi
export MACOSX_DEPLOYMENT_TARGET=11.0

export ENVROOT
export ROOUNFOLD_OLD
export CONDA_PREFIX="${ENVROOT}"
export ROOT_CXXMODULES="${ROOT_CXXMODULES:-ON}"

# Match the fish helper: if the env is in a split RooUnfold state, repair it by
# symlinking headers and libraries to the known-good external install.
if [[ ! -L "${ENVROOT}/include/RooUnfold.h" ]]; then
  mkdir -p "${ENVROOT}/include/RooUnfold_env_backup"
  mkdir -p "${ENVROOT}/lib/RooUnfold_env_backup"

  shopt -s nullglob
  for f in "${ENVROOT}"/include/RooUnfold*.h "${ENVROOT}"/include/RooUnfold*.tpp; do
    mv "${f}" "${ENVROOT}/include/RooUnfold_env_backup/"
  done

  for f in "${ROOUNFOLD_OLD}"/src/RooUnfold*.h; do
    ln -sf "${f}" "${ENVROOT}/include/$(basename "${f}")"
  done
  shopt -u nullglob

  [[ -e "${ROOUNFOLD_OLD}/libRooUnfold.so" ]] && \
    ln -sf "${ROOUNFOLD_OLD}/libRooUnfold.so" "${ENVROOT}/lib/libRooUnfold.so"
  [[ -e "${ROOUNFOLD_OLD}/libRooUnfold.rootmap" ]] && \
    ln -sf "${ROOUNFOLD_OLD}/libRooUnfold.rootmap" "${ENVROOT}/lib/libRooUnfold.rootmap"
  if [[ -e "${ROOUNFOLD_OLD}/RooUnfoldDict_rdict.pcm" ]]; then
    ln -sf "${ROOUNFOLD_OLD}/RooUnfoldDict_rdict.pcm" "${ENVROOT}/lib/RooUnfoldDict_rdict.pcm"
    ln -sf "${ROOUNFOLD_OLD}/RooUnfoldDict_rdict.pcm" "${ENVROOT}/lib/libRooUnfold_rdict.pcm"
  fi
fi

unset ROOT_INCLUDE_PATH
unset CPATH
unset DYLD_LIBRARY_PATH
unset LD_LIBRARY_PATH
unset ROOT_MODULE_PATH

export DYLD_FALLBACK_LIBRARY_PATH="${ROOUNFOLD_OLD}:${CONDA_PREFIX}/lib:/usr/local/lib:/usr/lib"
export DYLD_LIBRARY_PATH="${ROOUNFOLD_OLD}:${CONDA_PREFIX}/lib"
export LD_LIBRARY_PATH="${ROOUNFOLD_OLD}:${CONDA_PREFIX}/lib"
export LIBRARY_PATH="${ROOUNFOLD_OLD}:${CONDA_PREFIX}/lib"
export ROOT_INCLUDE_PATH="${ROOUNFOLD_OLD}/src:${CONDA_PREFIX}/include"
export CPATH="${ROOUNFOLD_OLD}/src:${CONDA_PREFIX}/include"
export ROOT_MODULE_PATH="${ROOUNFOLD_OLD}"
export PATH="${ENVROOT}/bin:/usr/bin:/bin:/usr/sbin:/sbin:${PATH:-}"

if [[ $# -eq 0 ]]; then
  exec "${ENVROOT}/bin/root"
else
  exec "$@"
fi
