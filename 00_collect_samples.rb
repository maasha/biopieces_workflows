#!/usr/bin/env ruby

# Take HiSeq pair-end files and merge forward reads into one file per sample,
# similarly for the reverse reads. Resulting files are loceted in the out_dir,
# prefixed with sample ID.

# Example of sampels.txt:
# 1       a/Sample_PRI_AARR_1/PRI_AARR_1_ATTACTCG-TATAGCCT_L001_R1_001.fastq.gz   a/Sample_PRI_AARR_1/PRI_AARR_1_ATTACTCG-TATAGCCT_L001_R2_001.fastq.gz
# 1       a/Sample_PRI_AARR_1/PRI_AARR_1_ATTACTCG-TATAGCCT_L001_R1_002.fastq.gz   a/Sample_PRI_AARR_1/PRI_AARR_1_ATTACTCG-TATAGCCT_L001_R2_002.fastq.gz
# 1       a/Sample_PRI_AARR_1/PRI_AARR_1_ATTACTCG-TATAGCCT_L002_R1_001.fastq.gz   a/Sample_PRI_AARR_1/PRI_AARR_1_ATTACTCG-TATAGCCT_L002_R2_001.fastq.gz
# 1       a/Sample_PRI_AARR_1/PRI_AARR_1_ATTACTCG-TATAGCCT_L002_R1_002.fastq.gz   a/Sample_PRI_AARR_1/PRI_AARR_1_ATTACTCG-TATAGCCT_L002_R2_002.fastq.gz
# 2       a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L001_R1_001.fastq.gz   a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L001_R2_001.fastq.gz
# 2       a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L001_R1_002.fastq.gz   a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L001_R2_002.fastq.gz
# 2       a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L001_R1_003.fastq.gz   a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L001_R2_003.fastq.gz
# 2       a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L002_R1_001.fastq.gz   a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L002_R2_001.fastq.gz
# 2       a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L002_R1_002.fastq.gz   a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L002_R2_002.fastq.gz
# 2       a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L002_R1_003.fastq.gz   a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L002_R2_003.fastq.gz
# 2       a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L002_R1_004.fastq.gz   a/Sample_PRI_AARR_2/PRI_AARR_2_TCCGGAGA-TATAGCCT_L002_R2_004.fastq.gz
# ...

require 'parallel'
require 'biopieces'
require 'csv'

cpus    = 20
samples = CSV.read("samples.txt", col_sep: "\s")
out_dir = "collected"
sids    = samples.each_with_object({}) { |e, a| a[e.first] = true }.keys

Parallel.each(sids, in_processes: cpus) do |sid|
  $stderr.puts "Collecting #{sid}"
  set = samples.select { |sample| sample.first == sid }

  r1_files = set.each_with_object([]) { |e, a| a << e[1] }
  r2_files = set.each_with_object([]) { |e, a| a << e[2] }

  BP.new.
  read_fastq(input: r1_files, encoding: :base_33).
  write_fastq(output: "#{sid}_R1.fastq").
  run

  BP.new.
  read_fastq(input: r2_files, encoding: :base_33).
  write_fastq(output: "#{sid}_R2.fastq").
  run

  $stderr.puts "Done"
end
