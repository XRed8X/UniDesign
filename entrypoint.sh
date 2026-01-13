#!/bin/bash
set -e

UNIDESIGN_BIN="/opt/unidesign/source/UniDesign"

if [ $# -eq 0 ]; then
    exec "$UNIDESIGN_BIN" --help
else
    exec "$UNIDESIGN_BIN" "$@"
fi