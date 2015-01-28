#!/usr/bin/env ruby

require 'biopieces'
require 'parallel'

cpus    = 10
out_dir = "Result-2015-01-28"

forward = "TGCCGTCTTCTGCTTG" # i7
reverse = "AATGATACGGCGACCA" # i5

pairs = []
Dir["Fastq/*Slaughterhouse*"].each_slice(2) { |pair| pairs << pair.sort }

Parallel.each(pairs, in_processes: cpus) do |pair|
  sample = File.basename(pair[0]).split('-').first

  $stderr.puts "Processing #{sample}"

  BP.new.
  read_fastq(input: pair[0], input2: pair[1], encoding: :base_33).
  plot_scores(output: "scores_preclip.png", terminal: :png, force: true).
  clip_primer(primer: forward, direction: :forward, mismatch_percent: 20).
  plot_histogram(key: :CLIP_PRIMER_POS, terminal: :png, output: "clip_primer_pos_forward.png", force: true).
  clip_primer(primer: reverse, reverse_complement: true, direction: :forward, mismatch_percent: 20).
  plot_histogram(key: :CLIP_PRIMER_POS, terminal: :png, output: "clip_primer_pos_reverse.png", force: true).
  trim_primer(primer: forward, direction: :forward, mismatch_percent: 20, overlap_min: 1).
  trim_primer(primer: reverse, reverse_complement: true, direction: :forward, mismatch_percent: 20, overlap_min: 1).
  plot_histogram(key: :SEQ_LEN, output: "lendist_postclip.png", terminal: :png, force: true).
  plot_scores(output: "scores_postclip.png", terminal: :png, force: true, count: true).
  trim_seq.
  merge_pair_seq.
  grab(evaluate: ":SEQ_LEN_LEFT >= 50").
  grab(evaluate: ":SEQ_LEN_RIGHT >= 50").
  mean_scores.
  grab(evaluate: ":SCORES_MEAN >= 25").
  mean_scores(local: true).
  grab(evaluate: ":SCORES_MEAN_LOCAL >= 15").
  split_pair_seq.
  random(number: 500_000, pairs: true).
  plot_scores(output: "scores_posttrim.png", terminal: :png, force: true, count: true).
  plot_histogram(key: :SEQ_LEN, output: "lendist_posttrim.png", terminal: :png, force: true).
  assemble_seq_spades(kmers: [21, 33, 55, 77, 99, 127], careful: true).
  sort(key: :SEQ_LEN).
  write_table(header: true, keys: [:SEQ_NAME, :SEQ_LEN], output: "contig_len.csv", force: true, pretty: true, commify: true).
  write_fasta(output: "contigs.fna", force: true).
  run(progress: false, output_dir: File.join(out_dir, sample), report: "assembly.html", email: "mail@maasha.dk")
end
