#!/bin/env zsh

# math through awk #{{{
calc() {
  # $1 is expression to calculate
  # $2 is number of decimal places (optional). Defaults to one if none given
  local df=${2:-1}
  awk "BEGIN { printf \"%.*f\n\", ${df}, ($1) }"
} #}}}
