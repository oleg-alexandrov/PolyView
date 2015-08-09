#!/bin/sh

self=$$
trap 'exit 1' TERM

TOPLEVEL="$(cd $(dirname $0)/.. && pwd)"
LIBEXEC="${TOPLEVEL}/libexec"

PROGRAM="${LIBEXEC}/$(basename $0)"

    case $(uname -s) in
        Linux)
        export LD_LIBRARY_PATH="$TOPLEVEL/lib"
        ;;
        Darwin)
        export DYLD_FALLBACK_LIBRARY_PATH="$TOPLEVEL/lib"
        export DYLD_FALLBACK_FRAMEWORK_PATH="$TOPLEVEL/lib"
        ;;
        *)
        die "Unknown OS: $(uname -s)"
        ;;
    esac

if ! test -f "${PROGRAM}"; then
    echo "Could not find ${PROGRAM}"
    exit 1
fi
exec "${PROGRAM}" "$@"
