#!/usr/bin/env ruby

require 'biopieces'
require 'csv'
require 'parallel'

cpus     = 20
forward  = "CCTAYGGGRBGCASCAG"     # 341
#forward  = "GTGCCAGCMGCCGCGGTAA"   # 515
reverse  = "GGACTACHVGGGTWTCTAAT"  # 806
out_dir  = "Amplicon-2015-09-07"

run_name = File.expand_path(__FILE__).split(File::SEPARATOR)[-2]
samples  = CSV.read("samples.txt", col_sep: "\s")

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start cleaning and dereplicating #{sample[0]}"

  p3 = BP.new.
  read_fastq(input: sample[1], input2: sample[2], encoding: :base_33).
  merge_pair_seq.
  plot_scores(terminal: :png, count: true, output: "p3_scores_pretrim_#{sample[0]}.png", force: true).
  split_pair_seq.
  clip_primer(primer: forward,  direction: :forward, mismatch_percent: 20, search_distance: 50).
  plot_histogram(key: :CLIP_PRIMER_POS, terminal: :png, output: "p3_clip_primer_pos_forward_#{sample[0]}.png", force: true).
  clip_primer(primer: reverse, direction: :forward, mismatch_percent: 20, search_distance: 100).
  plot_histogram(key: :CLIP_PRIMER_POS, terminal: :png, output: "p3_clip_primer_pos_reverse_#{sample[0]}.png", force: true).
  trim_primer(primer: forward,  direction: :forward, mismatch_percent: 20, overlap_min: 1).
  trim_primer(primer: reverse, direction: :forward, mismatch_percent: 20, overlap_min: 1).
  trim_seq.
  plot_histogram(key: :SEQ_LEN, terminal: :png, output: "p3_lendist_posttrim_#{sample[0]}.png", force: true).
  plot_scores(terminal: :png, count: true, output: "p3_scores_posttrim_#{sample[0]}.png", force: true).
  assemble_pairs(overlap_min: 1, mismatch_percent: 40, merge_unassembled: true, reverse_complement: true).
  grab(evaluate: ":SEQ_LEN >= 100").
  plot_histogram(key: :OVERLAP_LEN, terminal: :png, output: "p3_overlap_len_#{sample[0]}.png", force: true).
  plot_histogram(key: :HAMMING_DIST, terminal: :png, output: "p3_hamming_dist_#{sample[0]}.png", force: true).
  plot_residue_distribution(terminal: :png, count: true, output: "p3_residue_dist_postassembly_#{sample[0]}.png", force: true).
  mean_scores.
  grab(evaluate: ":SCORES_MEAN >= 25").
  mean_scores(local: true).
  grab(evaluate: ":SCORES_MEAN_LOCAL >= 15").
  plot_histogram(key: :SEQ_LEN, terminal: :png, output: "p3_lendist_postassembly_#{sample[0]}.png", force: true).
  plot_scores(terminal: :png, count: true, output: "p3_scores_postassembly_#{sample[0]}.png", force: true).
  write_fasta(output: "p3_clean_#{sample[0]}.fna", force: true).
  dereplicate_seq.
  write_table(output: "p3_derep_#{sample[0]}.tab", force: true, header: true, keys: [:SEQ_NAME, :SEQ, :SEQ_COUNT])

  p3.run(progress: false, verbose: false, output_dir: out_dir, report: "p3_#{sample[0]}.html")

  $stderr.puts "Done cleaning and dereplicating #{sample[0]}"
end

$stderr.puts "Start OTU clustering and chimera filtering"

p4 = BP.new.
read_table(input: "#{out_dir}/p3_derep*.tab", delimiter: "\t").
grab(evaluate: ":SEQ_COUNT > 1").
sort(key: :SEQ_COUNT, reverse: true).
cluster_otus(identity: 0.97).
uchime_ref.
add_key(key: :SEQ_NAME, prefix: "OTU_").
write_fasta(output: "p4_otus.fna", force: true).
classify_seq(prefix: '341F_806R_NR99').
grab(exact: true, keys: :RECORD_TYPE, select: "taxonomy").
write_table(output: "p4_classification_table.txt", header: true, force: true, keys: [:SEQ_NAME, :TAXONOMY])

p4.run(progress: true, verbose: false, output_dir: out_dir, report: "p4.html")

$stderr.puts "Done OTU clustering and chimera filtering"

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start remapping samples to OTUs for #{sample[0]}"

  p5 = BP.new.
  read_table(input: "#{out_dir}/p3_derep_#{sample[0]}.tab", delimiter: "\t").
  merge_values(keys: [:SEQ_NAME, :SEQ_COUNT], delimiter: ":count=").
  usearch_global(database: "#{out_dir}/p4_otus.fna", identity: 0.97, strand: "plus").
  grab(exact: true, keys: :TYPE, select: 'H').
  add_key(key: :SAMPLE, value: sample[0]).
  write_table(output: "p5_usearch_global_#{sample[0]}.tab", header: true, force: true, keys: [:TYPE, :Q_ID, :S_ID, :SAMPLE])

  p5.run(progress: false, verbose: false, output_dir: out_dir, report: "p5_#{sample[0]}.html")

  $stderr.puts "Done remapping samples to OTUs for #{sample[0]}"
end

$stderr.puts "Start collecting OTU table"

p6 = BP.new.
read_table(input: "#{out_dir}/p5_usearch_global_*.tab", delimiter: "\t").
split_values(key: :Q_ID, keys: [:Q_ID, :SEQ_COUNT], delimiter: ":count=").
collect_otus.
grab(exact: true, keys: :RECORD_TYPE, select: 'OTU').
merge_table(input: "#{out_dir}/p4_classification_table.txt", key: :OTU, keys: [:OTU, :TAXONOMY], delimiter: "\t").
collapse_otus.
plot_heatmap(skip: [:RECORD_TYPE, :OTU, :TAXONOMY],terminal: :png, output: "p6_heatmap.png", force: true, xlabel: "Samples", ylabel: "OTUs", logscale: true).
write_table(header: true, output: "p6_otu_table.txt", skip: [:RECORD_TYPE], force: true)

p6.run(progress: true, verbose: false, output_dir: out_dir, report: "p6.html")

$stderr.puts "Done collecting OTU table"

$stderr.puts "Start creating tree"

p7 = BP.new.
read_fasta(input: "#{out_dir}/p4_otus.fna").
align_seq_mothur.
degap_seq(columns_only: true).
write_fasta(output: "p7_aligned.fna", force: true).
write_tree(output: "p7.tree", force: true).
run(verbose: true, output_dir: out_dir, report: "p7.html", email: "mail@maasha.dk", subject: "#{run_name} done")

$stderr.puts "Done creating tree"
