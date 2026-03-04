#!/usr/bin/env python

import argparse
import json

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
  
def main(batch_id, element, json_path):
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
  parser = argparse.ArgumentParser()
  parser.add_argument('--batch_id', default='', help='Task id or batch number from which the element are retrieved.')
  parser.add_argument(
    '--element',
    default='all',
    choices=['all', 'start', 'end', 'model'],
    help='Select the element (start, end or model) that is retrieved by the script.')
  parser.add_argument('--json_path', default='./batches.json', help='Path of the json from which the element are retrieved.')

  parsed = parser.parse_args()
  main(
    parsed.batch_id,
    parsed.element,
    parsed.json_path
  )
