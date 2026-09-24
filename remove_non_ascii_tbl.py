# remove_non_ascii_tbl.py

import sys
import unicodedata

def remove_non_ascii(text):
    """
    Convert non-ASCII characters to closest ASCII equivalent.
    If no equivalent, remove the character.
    """
    normalized_text = unicodedata.normalize('NFKD', text)
    ascii_text = normalized_text.encode('ascii', 'ignore').decode('ascii')
    return ascii_text

def process_tbl_file(input_file, output_file):
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    cleaned_lines = [remove_non_ascii(line) for line in lines]

    with open(output_file, 'w', encoding='ascii') as f:
        f.writelines(cleaned_lines)

    print(f"Processed .tbl file saved as: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python remove_non_ascii_tbl.py input_file.tbl output_file.tbl")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_tbl_file(input_file, output_file)
