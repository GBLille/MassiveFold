#!/usr/bin/env python3
"""CLI skeleton for MassiveFold."""

import argparse
import sys

def _build_parser():
  parser = argparse.ArgumentParser(prog="massivefold")
  subparsers = parser.add_subparsers(dest="command")

  subparsers.add_parser("run", help="Run MassiveFold (skeleton)")
  subparsers.add_parser("screening", help="Run MassiveFold screening (skeleton)")

  return parser


def _print_skeleton_message(command, forwarded_args):
  print(f"massivefold {command}: CLI skeleton only (no runtime implementation yet).")
  if forwarded_args:
    print("Forwarded arguments:", " ".join(forwarded_args))
  return 0

def main(argv=None):
  parser = _build_parser()
  args, unknown = parser.parse_known_args(argv)

  if not args.command:
    parser.print_help()
    return 0

  if args.command == "run":
    return _print_skeleton_message("run", unknown)
  if args.command == "screening":
    return _print_skeleton_message("screening", unknown)

  parser.print_help()
  return 1

if __name__ == "__main__":
  raise SystemExit(main(sys.argv[1:]))
