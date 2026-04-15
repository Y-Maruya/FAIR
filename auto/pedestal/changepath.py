#!/usr/bin/env python3
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description='Replace placeholders in YAML configuration file'
    )
    parser.add_argument('yaml_file', help='Path to YAML file')
    parser.add_argument('-i', '--input', help='Input file path')
    parser.add_argument('-o', '--output', help='Output directory path')
    parser.add_argument('-r', '--run-number', help='Run number')
    
    args = parser.parse_args()
    
    try:
        # Read the YAML file
        with open(args.yaml_file, 'r') as f:
            content = f.read()
        
        # Replace placeholders
        if args.input:
            content = content.replace('###input-file###', args.input)
        if args.output:
            content = content.replace('###out-file###', args.output)
        if args.run_number:
            content = content.replace('###run-number###', args.run_number)
        
        # Write back to the YAML file
        with open(args.yaml_file, 'w') as f:
            f.write(content)
        
        print(f"Successfully updated {args.yaml_file}")
        
    except FileNotFoundError:
        print(f"Error: File '{args.yaml_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()