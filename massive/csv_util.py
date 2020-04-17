import csv
import re

#Headers: Removes everything inside the parenthesis and stripes all spaces.
def parse_csv(csv_file_path, delimiter=','):
    try:
        csv_file = open(csv_file_path, 'rU', encoding='utf-8-sig')
        elements = get_csv_elements(csv_file, delimiter)
    except UnicodeDecodeError:
        csv_file = open(csv_file_path, 'rU', encoding='latin-1')
        elements = get_csv_elements(csv_file, delimiter)

    csv_file.close()
    return elements

def get_csv_elements(file_descriptor, delimiter):
    delimiter_dict = {";" : ",", "," : ";"}
    csv_reader = csv.DictReader(file_descriptor, delimiter=delimiter)
    if len(csv_reader.fieldnames) < 2:
        delimiter=delimiter_dict[delimiter]
        file_descriptor.seek(0)
        csv_reader = csv.DictReader(file_descriptor, delimiter=delimiter)
    return [{re.sub(r'(\s+|\([^)]*\))', '', k):v for k, v in row.items()} for row in csv_reader], [re.sub(r'(\s+|\([^)]*\))', '', fieldname) for fieldname in csv_reader.fieldnames]

def write_csv(csv_file_path, fieldnames, dict_list):
    with open(csv_file_path, 'w', encoding='utf-8-sig') as csv_file:
        csv_writer = csv.DictWriter(csv_file, fieldnames, delimiter=';')
        csv_writer.writeheader()
        for d in dict_list:
            csv_writer.writerow(d)
    return csv_file_path
