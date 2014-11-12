#!/usr/bin/env ruby

require 'biopieces'
require 'csv'
require 'parallel'

cpus    = 4
forward = 'GCGGCTTATCARRTNGARGG'
reverse = 'GTGTAACACYMNGGYTCRTT'
out_dir = 'fungi_2014-10-22'

samples = CSV.read("fungi.txt", col_sep: "\s")

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start cleaning and dereplicating #{sample[0]}"

  p1 = BP.new.
  read_fastq(input: sample[1], input2: sample[2], encoding: :base_33).
  plot_scores(terminal: :png, count: true, output: "scores_pretrim_#{sample[0]}.png", force: true).
  clip_primer(primer: forward,  direction: :forward, mismatch_percent: 20, search_distance: 50).
  clip_primer(primer: forward,  direction: :reverse, reverse_complement: true, mismatch_percent: 20, search_distance: 50).
  plot_histogram(key: :CLIP_PRIMER_POS, terminal: :png, output: "clip_primer_pos_forward_#{sample[0]}.png", force: true).
  clip_primer(primer: reverse, direction: :forward, mismatch_percent: 20, search_distance: 50).
  clip_primer(primer: reverse, direction: :reverse, reverse_complement: true, mismatch_percent: 20, search_distance: 50).
  plot_histogram(key: :CLIP_PRIMER_POS, terminal: :png, output: "clip_primer_pos_reverse_#{sample[0]}.png", force: true).
  trim_primer(primer: forward, direction: :forward, mismatch_percent: 20, overlap_min: 1).
  trim_primer(primer: forward, direction: :reverse, reverse_complement: true, mismatch_percent: 20, overlap_min: 1).
  trim_primer(primer: reverse, direction: :forward, mismatch_percent: 20, overlap_min: 1).
  trim_primer(primer: reverse, direction: :reverse, reverse_complement: true, mismatch_percent: 20, overlap_min: 1).
  trim_seq.
  plot_histogram(key: :SEQ_LEN, terminal: :png, output: "lendist_posttrim_#{sample[0]}.png", force: true).
  plot_scores(terminal: :png, count: true, output: "scores_posttrim_#{sample[0]}.png", force: true).
  assemble_pairs(overlap_min: 40, mismatch_percent: 40, reverse_complement: true).
  mean_scores.
  grab(evaluate: ":SCORES_MEAN >= 25").
  mean_scores(local: true).
  grab(evaluate: ":SCORES_MEAN_LOCAL >= 15").
  plot_histogram(key: :SEQ_LEN, terminal: :png, output: "lendist_postassembly_#{sample[0]}.png", force: true).
  plot_scores(terminal: :png, count: true, output: "scores_postassembly_#{sample[0]}.png", force: true).
  write_fasta(output: "clean_#{sample[0]}.fna", force: true).
  dereplicate_seq.
  write_table(output: "derep_#{sample[0]}.tab", force: true, header: true, keys: [:SEQ_NAME, :SEQ, :SEQ_COUNT])

  p1.run(progress: false, verbose: false, output_dir: out_dir, report: "p1_#{sample[0]}.html")

  $stderr.puts "Done cleaning and dereplicating #{sample[0]}"
end

$stderr.puts "Start OTU clustering and chimera filtering"

p2 = BP.new.
read_table(input: "#{out_dir}/derep*.tab", delimiter: "\t").
grab(evaluate: ":SEQ_COUNT > 1").
sort(key: :SEQ_COUNT, reverse: true).
cluster_otus(identity: 0.85).
add_key(key: :SEQ_NAME, prefix: "OTU_").
write_fasta(output: "otus.fna", force: true)

p2.run(progress: true, verbose: false, output_dir: out_dir, report: "p2.html")

$stderr.puts "Done OTU clustering and chimera filtering"

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start remapping samples to OTUs for #{sample[0]}"

  p3 = BP.new.
  read_table(input: "#{out_dir}/derep_#{sample[0]}.tab", delimiter: "\t").
  merge_values(keys: [:SEQ_NAME, :SEQ_COUNT], delimiter: ":count=").
  usearch_global(database: "#{out_dir}/otus.fna", identity: 0.85, strand: "plus").
  grab(exact: true, keys: :TYPE, select: 'H').
  add_key(key: :SAMPLE, value: sample[0]).
  write_table(output: "usearch_global_#{sample[0]}.tab", header: true, force: true, keys: [:TYPE, :Q_ID, :S_ID, :SAMPLE])

  p3.run(progress: false, verbose: false, output_dir: out_dir, report: "p3_#{sample[0]}.html")

  $stderr.puts "Done remapping samples to OTUs for #{sample[0]}"
end

$stderr.puts "Start collecting OTU table"

p4 = BP.new.
read_table(input: "#{out_dir}/usearch_global_*.tab", delimiter: "\t").
split_values(key: :Q_ID, keys: [:Q_ID, :SEQ_COUNT], delimiter: ":count=").
collect_otus.
grab(exact: true, keys: :RECORD_TYPE, select: 'OTU').
write_table(header: true, output: "otu_table.tsv", skip: [:RECORD_TYPE])

p4.run(progress: true, verbose: false, output_dir: out_dir, report: "p4.html")

$stderr.puts "Done collecting OTU table"
