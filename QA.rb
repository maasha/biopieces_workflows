#!/usr/bin/env ruby

require 'biopieces'
require 'csv'
require 'parallel'

cpus     = 20
out_dir  = "QA_2015-07-10"
run_name = File.expand_path(__FILE__).split(File::SEPARATOR)[-2]
samples  = CSV.read("samples.txt", col_sep: "\s")

Parallel.each(samples, in_processes: cpus) do |sample|
  $stderr.puts "Start QA: #{sample[0]}"

  p1 = BP.new.
  read_fastq(input: sample[1], input2: sample[2], encoding: :base_33).
  merge_pair_seq.
  plot_scores(terminal: :png, count: true, output: "p1_scores_#{sample[0]}.png", force: true).
  split_pair_seq.
  count.
  add_key(key: :SAMPLE, value: sample[0]).
  grab(exact: true, keys: :RECORD_TYPE, select: 'count').
  write_table(header: true, output: "p1_#{sample[0]}_count.tab", force: true)

  p1.run(progress: false, output_dir: out_dir, report: "p1_#{sample[0]}.html")

  $stderr.puts "Done QA: #{sample[0]}"
end

$stderr.puts "Start collecting QA"

p2 = BP.new.
read_table(input: "#{out_dir}/p1*_count.tab", delimiter: "\t").
sort(key: :COUNT, reverse: true).
plot_histogram(key: :SAMPLE, value: :COUNT, output: "p2_count.png", terminal: :png, force: true).
write_table(header: true, pretty: true, commify: true, output: "p2_count.tab", force: true, skip: [:RECORD_TYPE])

p2.run(progress: true, output_dir: out_dir, report: "p2.html")

$stderr.puts "Done collecting QA"
