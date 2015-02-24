#!/usr/bin/env ruby

require 'biopieces'
require 'parallel'

cpus    = 10
out_dir = "Result-2015-02-24"

# GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG # i7
# TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG  # i5
forward = "GATGTGTATAAGAGACAG" # i7
reverse = "GATGTGTATAAGAGACAG" # i5

pairs = []
Dir["Fastq/*Slaughterhouse*"].each_slice(2) { |pair| pairs << pair.sort }

Parallel.each(pairs, in_processes: 20) do |pair|
  sample = File.basename(pair[0]).split('-').first

  $stderr.puts "Processing #{sample}"

  p1 = BP.new.
  read_fastq(input: pair[0], input2: pair[1], encoding: :base_33).
  plot_scores(output: "scores_preclip.png", terminal: :png, force: true).
  plot_nucleotide_distribution(output: "nucdist_preclip.png", terminal: :png, force: true).
  clip_primer(primer: forward, reverse_complement: false, direction: :forward, mismatch_percent: 20).
  plot_histogram(key: :CLIP_PRIMER_POS, terminal: :png, output: "clip_primer_pos_forward.png", force: true).
  clip_primer(primer: reverse, reverse_complement: true, direction: :forward, mismatch_percent: 20).
  plot_histogram(key: :CLIP_PRIMER_POS, terminal: :png, output: "clip_primer_pos_reverse.png", force: true).
  trim_primer(primer: forward, reverse_complement: false, direction: :forward, mismatch_percent: 20, overlap_min: 1).
  trim_primer(primer: reverse, reverse_complement: true, direction: :forward, mismatch_percent: 20, overlap_min: 1).
  plot_histogram(key: :SEQ_LEN, output: "lendist_postclip.png", terminal: :png, force: true).
  plot_scores(output: "scores_postclip.png", terminal: :png, force: true, count: true).
  plot_nucleotide_distribution(output: "nucdist_postclip.png", terminal: :png, force: true, count: true).
  trim_seq.
  merge_pair_seq.
  grab(evaluate: ":SEQ_LEN_LEFT >= 50").
  grab(evaluate: ":SEQ_LEN_RIGHT >= 50").
  mean_scores.
  grab(evaluate: ":SCORES_MEAN >= 25").
  mean_scores(local: true, window_size: 10).
  grab(evaluate: ":SCORES_MEAN_LOCAL >= 10").
  usearch_global(database: "/home/maasha/data/UniVec_Core.txt", identity: 0.9, strand: "both").
  split_pair_seq.
  plot_scores(output: "scores_posttrim.png", terminal: :png, force: true, count: true).
  plot_histogram(key: :SEQ_LEN, output: "lendist_posttrim.png", terminal: :png, force: true).
  plot_nucleotide_distribution(output: "nucdist_postclean.png", terminal: :png, force: true, count: true).
  trim_seq.
  write_fasta(output: "reads_clean.fna", force: true).
  grab(select: "S_ID", keys_only: true).
  write_table(output: "contaminants.txt", force: true, header: true).
  write_table(output: "IDs_to_remove.txt", keys: [:Q_ID], force: true, header: true)

  p1.run(progress: false, output_dir: File.join(out_dir, sample), report: "cleaning.html")

  p2 = BP.new.
  read_fasta(input: File.join(out_dir, sample, "reads_clean.fna")).
  merge_pair_seq.
  grab(reject_file: File.join(out_dir, sample, "IDs_to_remove.txt"), exact: true, keys: [:SEQ_NAME]).
  split_pair_seq.
  assemble_seq_spades(kmers: [21, 33, 55, 77, 99, 127], careful: true).
  grab(evaluate: ":SEQ_LEN >= 1000").
  sort(key: :SEQ_LEN).
  write_table(header: true, keys: [:SEQ_NAME, :SEQ_LEN], output: "contig_len.txt", force: true, pretty: true, commify: true).
  write_fasta(output: "contigs.fna", force: true)

  p2.run(progress: false, output_dir: File.join(out_dir, sample), report: "assembly.html", email: "mail@maasha.dk", subject: "#{sample} done")

  $stderr.puts "Processing #{sample} done."

  exit
end

