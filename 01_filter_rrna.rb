#!/usr/bin/env ruby

# Identify all reads NOT matching ribosomal RNAs and write IDs to file.
# Next the IDs are uniqued.
# Read in all reads in interleaved order and select all mathing to the list of
# unique IDs and save these reads as non-ribosomal.
# Read in all reads in interleaved order and reject all mathing to the list of
# unique IDs and save these reads as ribosomal.

# Example of samples_collected.txt:
# 1	collected/1_R1.fastq	collected/1_R2.fastq
# 2	collected/2_R1.fastq	collected/2_R2.fastq
# 3	collected/3_R1.fastq	collected/3_R2.fastq
# 4	collected/4_R1.fastq	collected/4_R2.fastq
# 5	collected/5_R1.fastq	collected/5_R2.fastq
# 6	collected/6_R1.fastq	collected/6_R2.fastq
# 7	collected/7_R1.fastq	collected/7_R2.fastq
# 9	collected/9_R1.fastq	collected/9_R2.fastq
# 10	collected/10_R1.fastq	collected/10_R2.fastq
# ...

require 'parallel'
require 'biopieces'
require 'csv'

cpus    = 20
samples = CSV.read("samples_collected.txt", col_sep: "\s")
files   = samples.each_with_object([]) { |e, a| a << e[1] && a << e[2] }
out_dir = "filter_rrna"

Parallel.each(files, in_processes: cpus) do |file|
  basename = File.basename(file, ".fastq")
  $stderr.puts "Start locating rRNAs #{basename}"

  BP.new.
  read_fastq(encoding: :base_33, input: file).
  filter_rrna.
  write_table(output: "#{basename}_ids.tab", keys: [:SEQ_NAME]).
  run(output_dir: out_dir, report: "#{basename}_sortmerna.html")

  $stderr.puts "Done locating rRNAs #{basename}"
end

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start compiling filter id list for #{sample[0]}"

  BP.new.
  read_table(input: [File.join(out_dir, "#{sample[0]}_R1_ids.tab"), File.join(out_dir, "#{sample[0]}_R2_ids.tab")], delimiter: "\t").
  split_values(key: :V0, delimiter: ' ').
  unique_values(key: :V0_0).
  write_table(output: "#{sample[0]}_ids.tab", keys: [:V0_0]).
  run(output_dir: out_dir, report: "#{sample[0]}_filter_list.html")

  $stderr.puts "Done compiling filter id list for #{sample[0]}"
end

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start isolating non-rRNA reads #{sample[0]}"

  BP.new.
  read_fastq(encoding: :base_33, input: sample[1], input2: sample[2]).
  merge_pair_seq.
  split_values(key: :SEQ_NAME, delimiter: ' ').
  grab(select_file: File.join(out_dir, "#{sample[0]}_ids.tab"), keys: :SEQ_NAME_0, exact: true).
  split_pair_seq.
  write_fastq(output: "#{sample[0]}_norrna.fq").
  run(output_dir: out_dir, report: "#{sample[0]}_norrna.html")

  $stderr.puts "Done isolating non-rRNA reads #{sample[0]}"
end

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start isolating rRNA reads #{sample[0]}"

  BP.new.
  read_fastq(encoding: :base_33, input: sample[1], input2: sample[2]).
  merge_pair_seq.
  split_values(key: :SEQ_NAME, delimiter: ' ').
  grab(reject_file: File.join(out_dir, "#{sample[0]}_ids.tab"), keys: :SEQ_NAME_0, exact: true).
  split_pair_seq.
  write_fastq(output: "#{sample[0]}_rrna.fq").
  run(output_dir: out_dir, report: "#{sample[0]}_rrna.html")

  $stderr.puts "Done isolating rRNA reads #{sample[0]}"
end
