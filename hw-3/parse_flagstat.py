#!/usr/bin/env python3
import sys
import re

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 parse_flagstat.py <flagstat.txt>")
        sys.exit(1)

    filename = sys.argv[1]
    
    with open(filename, 'r') as f:
        content = f.read()

    match = re.search(r'mapped \(([\d\.]+)%', content)
    if match:
        mapped_percent = float(match.group(1))
        print(mapped_percent)
    else:
        print("Error: Could not find mapped percentage in the output.")
        sys.exit(1)

if __name__ == "__main__":
    main()
