#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: \$0 <exec> <lcc_in> <sd_in> <results_path>"
    exit 1
fi

EXEC_PATH="$1"
LCC_PATH="$2"
SD_PATH="$3"
RESULTS_PATH="$4"

EXEC_OUTPUT=$("$EXEC_PATH" "$LCC_PATH" "$SD_PATH")
LAST_LINE=$(echo "$EXEC_OUTPUT" | tail -n 1)

echo $LAST_LINE

# Print the required information
{
    echo -e "lcc:\n[$(cat "$LCC_PATH")]\n\n"
    echo -e "sd:\n[$(cat "$SD_PATH")]\n\n"
    echo -e "$LAST_LINE\n-----------------------------\n"
} >> $RESULTS_PATH
