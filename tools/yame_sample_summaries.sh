#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: yame_sample_summaries.sh --cg FILE --aggregate-output FILE
       [--chromhmm-cm FILE --chromhmm-output FILE] [--yame COMMAND]

Aggregate an indexed per-pixel YAME matrix once and fan the stream out to the
sample-level CpG summary and, when requested, the ChromHMM summary.
EOF
}

yame=yame
cg=
aggregate_output=
chromhmm_cm=
chromhmm_output=
while (($#)); do
    case "$1" in
        --yame) yame=$2; shift 2 ;;
        --cg) cg=$2; shift 2 ;;
        --aggregate-output) aggregate_output=$2; shift 2 ;;
        --chromhmm-cm) chromhmm_cm=$2; shift 2 ;;
        --chromhmm-output) chromhmm_output=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done

if [[ -z $cg || -z $aggregate_output ]]; then
    usage
    exit 2
fi
if [[ -n $chromhmm_cm && -z $chromhmm_output ]]; then
    echo "--chromhmm-output is required with --chromhmm-cm" >&2
    exit 2
fi

mkdir -p "$(dirname "$aggregate_output")"
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/yame-summary.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT
mkfifo "$tmpdir/aggregate.fifo"

"$yame" summary - < "$tmpdir/aggregate.fifo" > "$tmpdir/aggregate.tsv" &
aggregate_pid=$!

if [[ -n $chromhmm_cm ]]; then
    if [[ ! -r $chromhmm_cm ]]; then
        echo "ChromHMM mask is not readable: $chromhmm_cm" >&2
        exit 1
    fi
    mkdir -p "$(dirname "$chromhmm_output")"
    mkfifo "$tmpdir/chromhmm.fifo"
    "$yame" summary -m "$chromhmm_cm" - \
        < "$tmpdir/chromhmm.fifo" > "$tmpdir/chromhmm.tsv" &
    chromhmm_pid=$!
fi

set +e
if [[ -n $chromhmm_cm ]]; then
    "$yame" rowop -o musum "$cg" \
        | tee "$tmpdir/aggregate.fifo" > "$tmpdir/chromhmm.fifo"
    stream_status=$?
else
    "$yame" rowop -o musum "$cg" > "$tmpdir/aggregate.fifo"
    stream_status=$?
fi
wait "$aggregate_pid"
aggregate_status=$?
chromhmm_status=0
if [[ -n $chromhmm_cm ]]; then
    wait "$chromhmm_pid"
    chromhmm_status=$?
fi
set -e

if ((stream_status || aggregate_status || chromhmm_status)); then
    echo "YAME summary failed: stream=$stream_status aggregate=$aggregate_status chromhmm=$chromhmm_status" >&2
    exit 1
fi

mv "$tmpdir/aggregate.tsv" "$aggregate_output"
if [[ -n $chromhmm_cm ]]; then
    mv "$tmpdir/chromhmm.tsv" "$chromhmm_output"
fi
