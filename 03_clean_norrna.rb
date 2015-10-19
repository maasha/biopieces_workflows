#!/usr/bin/env ruby

# Clip and trim primers from all non-ribosomal reads, trim according to
# qualilty scores, and plot some stats.

require 'biopieces'
require 'parallel'

cpus    = 20
out_dir = "clean_norrna"
samples = CSV.read("samples_collected.txt", col_sep: "\s")
forward = "TGCCGTCTTCTGCTTG" # i7
reverse = "AATGATACGGCGACCA" # i5

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start cleaning #{sample[0]}"

  BP.new.
  read_fastq(input: File.join('filter_rrna', "#{sample[0]}_norrna.fq"), encoding: :base_33).
  clip_primer(primer: forward, reverse_complement: false, direction: :forward, mismatch_percent: 20).
  plot_histogram(key: :CLIP_PRIMER_POS, terminal: :png, output: "#{sample[0]}_clip_forward.png", force: true).
  clip_primer(primer: reverse, reverse_complement: true, direction: :forward, mismatch_percent: 20).
  plot_histogram(key: :CLIP_PRIMER_POS, terminal: :png, output: "#{sample[0]}_clip_reverse.png", force: true).
  trim_primer(primer: forward, direction: :forward, mismatch_percent: 20, overlap_min: 1).
  trim_primer(primer: reverse, direction: :forward, mismatch_percent: 20, overlap_min: 1).
  trim_seq.
  merge_pair_seq.
  grab(evaluate: ":SEQ_LEN_LEFT >= 50").
  grab(evaluate: ":SEQ_LEN_RIGHT >= 50").
  split_pair_seq.
  plot_histogram(key: :SEQ_LEN, terminal: :png, output: "#{sample[0]}_lendist.png", force: true).
  plot_scores(terminal: :png, count: true, output: "#{sample[0]}_scores.png", force: true).
  write_fastq(output: "#{sample[0]}_clean.fq", force: true).
  run(output_dir: out_dir, report: "#{sample[0]}_clean.html")

  $stderr.puts "Done cleaning #{sample[0]}"
end
