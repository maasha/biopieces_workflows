#!/usr/bin/env ruby

# Read in cleaned non-ribosomal RNA matching reads and assemble read pairs
# allowing for unassembled reads and write a table of IDs for read pairs that
# were assembled.
#
# Next assmbled and singleton reads are separated.

require 'biopieces'
require 'parallel'

cpus    = 10
out_dir = "assemble_pairs_norrna"
samples = CSV.read("samples_collected.txt", col_sep: "\s")

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start assembling pairs for #{sample[0]}"

  BP.new.
  read_fastq(input: File.join("clean_norrna", "#{sample[0]}_clean.fq"), encoding: :base_33).
  assemble_pairs(overlap_min: 20, mismatch_percent: 40, allow_unassembled: true, reverse_complement: true).
  write_fasta(output: "#{sample[0]}_assemble_pairs_all.fna", force: true).
  grab(select: 'overlap', keys: [:SEQ_NAME]).
  plot_histogram(key: :OVERLAP_LEN, terminal: :png, output: "#{sample[0]}_overlap.png", force: true).
  plot_histogram(key: :HAMMING_DIST, terminal: :png, ylogscale: true, output: "#{sample[0]}_hamming.png", force: true).
  write_table(output: "#{sample[0]}_paired_ids.tab", keys: [:SEQ_NAME], force: true).
  run(progress: false, verbose: false, output_dir: out_dir, report: "#{sample[0]}_assemble_pairs.html")

  $stderr.puts "Done assembling pairs for #{sample[0]}"
end

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start Isolating paired reads for #{sample[0]}"

  BP.new.
  read_fasta(input: File.join(out_dir, "#{sample[0]}_assemble_pairs_all.fna")).
  grab(select_file: File.join(out_dir, "#{sample[0]}_paired_ids.tab"), keys: :SEQ_NAME, exact: true).
  plot_histogram(key: :SEQ_LEN, terminal: :png, output: "#{sample[0]}_pairs_lendist.png", force: true).
  write_fasta(output: "#{sample[0]}_norrna_assemble_pairs_pairs.fna", force: true).
  run(output_dir: out_dir, report: "#{sample[0]}_isolate_paired.html")

  $stderr.puts "Done isolating paired reads for #{sample[0]}"
end

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start Isolating singleton reads for #{sample[0]}"

  BP.new.
  read_fasta(input: File.join(out_dir, "#{sample[0]}_assemble_pairs_all.fna")).
  grab(reject_file: File.join(out_dir, "#{sample[0]}_paired_ids.tab"), keys: :SEQ_NAME, exact: true).
  plot_histogram(key: :SEQ_LEN, terminal: :png, output: "#{sample[0]}_singletons_lendist.png", force: true).
  write_fasta(output: "#{sample[0]}_norrna_assemble_pairs_singletons.fna", force: true).
  run(output_dir: out_dir, report: "#{sample[0]}_isolate_singletons.html")

  $stderr.puts "Done isolating singletons reads for #{sample[0]}"
end
