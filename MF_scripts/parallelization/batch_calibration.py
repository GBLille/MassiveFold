#!/usr/bin/env python

import os
from math import floor, ceil
from absl import flags, app

FLAGS = flags.FLAGS
flags.DEFINE_integer('wall_time', 20, 'Inference time in hour to not exceed.')
flags.DEFINE_string('logs_dir', '', 'Log directory path of the run.')
flags.DEFINE_float('add_excess', 0.1, 'Excess time proportion for the inference of a single prediction')

def extract_longer(jobarray_path):
  time_lines = os.popen(f"cat {jobarray_path}/jobarray_* | grep 'predict time'").read()
  i_thing = time_lines.split(' ')[0]
  time_lines_list = time_lines.split(i_thing)

  lines = [line[1:-1] for line in time_lines_list]
  times = [float(line.split(' ')[-1].replace('s', '')) for line in lines[1:]]
  return ceil(max(times))

def max_pred_nb_for_walltime(wall_time, single_pred_time):
  wall_time_m = wall_time*60
  wall_time_s = wall_time_m*60
  max_number = wall_time_s/single_pred_time
  return floor(max_number)

def safen_time(single_pred_time, excess_proportion):
  return floor(single_pred_time + (single_pred_time*excess_proportion))

def main(argv):
  path_to_logs = FLAGS.logs_dir
  wall_time_h = FLAGS.wall_time

  max_time = extract_longer(path_to_logs)
  safe_time = safen_time(max_time, FLAGS.add_excess)
  print(max_pred_nb_for_walltime(wall_time_h, safe_time))

if __name__ == "__main__":
  app.run(main)
