#!/usr/bin/env ruby

# Merge singleton non-ribosomal reads end-to-end.

require 'biopieces'
require 'parallel'

cpus    = 20
out_dir = "assemble_pairs_norrna"
samples = CSV.read("samples_collected.txt", col_sep: "\s")

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start merging pairs for #{sample[0]}"

  BP.new.
  read_fasta(input: File.join("assemble_pairs_norrna", "#{sample[0]}_norrna_assemble_pairs_singletons.fna")).
  merge_pair_seq.
  write_fasta(output: "#{sample[0]}_merge_singletons.fna", force: true).
  run(output_dir: out_dir)

  $stderr.puts "Done merging singletons reads for #{sample[0]}"
end
