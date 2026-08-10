#!/usr/bin/env zsh

# math through awk #{{{
calc() {
  # $1 is expression to calculate
  # $2 is number of decimal places (optional). Defaults to one if none given
  local df=${2:-1}
  awk "BEGIN { printf \"%.*f\n\", ${df}, ($1) }"
} #}}}
# zero-pad a number to a fixed width #{{{
zpad() {
  # $1 is the number to pad
  # $2 is the total width (optional). Defaults to four if none given
  # e.g., zpad 7 4 -> 0007; useful for numbered output files that should sort
  local width=${2:-4}
  printf "%0*d\n" "${width}" "$1"
} #}}}
