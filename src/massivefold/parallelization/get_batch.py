#!/usr/bin/env python

import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('--batch_id', default='', required=True, help='Task id or batch number from which the element are retrieved.')
parser.add_argument(
  '--element',
  default='all',
  choices=['all', 'start', 'end', 'model'],
  required=True,
  help='Select the element (start, end or model) that is retrieved by the script.')
parser.add_argument('--json_path', default='./batches.json', required=True, help='Path of the json from which the element are retrieved.')

def get_single_element(batch_id, element, json_path):
  with open(json_path, 'r') as f:
    all_batches = json.load(f)

  try:
    print(all_batches[batch_id][element])
  except KeyError:
    raise ValueError(f"Either no batch {batch_id} or no element '{element}' for this batch.")

def get_all_elements(all_batches, batch_id, element):
  for element in ['start', 'end', 'model']:
    try:
      print(all_batches[batch_id][element])
    except KeyError:
      raise ValueError(f"Either no batch {batch_id} or no element '{element}' for this batch.")

def main():
  args = parser.parse_args()
  batch_id = args.batch_id
  element = args.element
  json_path = args.json_path
  if not batch_id or not element:
    raise ValueError('\nUsage: ./get_batch.py --batch_id [str:batch_id] --element [str:element] [[options]]')
  if element not in ['start', 'end', 'model', 'all']:
    raise ValueError('--element is either start, end or model.')

  with open(json_path, 'r') as f:
    all_batches = json.load(f)

  if element != 'all':
    get_single_element(batch_id, element, json_path)
  else:
    get_all_elements(all_batches, batch_id, element)

if __name__ == "__main__":
  main()
