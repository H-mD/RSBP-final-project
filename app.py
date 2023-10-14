# from Bio import SeqIO
from flask import Flask, render_template
import json
import os
import time
import csv
import re
import uuid
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from collections import defaultdict
from decimal import Decimal

app = Flask(__name__)
genome_file = "./test_data/sequence.gb"

class CGViewBuilder:
    VERSION = '1.1.0'

    def __init__(self, sequence_path, options={}):
        self.map_id = options.get('map_id', str(uuid.uuid4()))
        self.map_name = options.get('map_name')
        self.cgview = self.initialize_cgview()
        self.options = options
        self.features = []
        self.plots = []
        self.tracks = []
        self.blast_tracks = []
        self.debug = False
        self.config = {}
        self.read_sequence(sequence_path)
        # self.build_genetic_code()
        self.common_genetic_code = None
        if options.get('analysis_path'):
            self.read_gff_analysis(options['analysis_path'])
        if options.get('blasts'):
            self.read_blasts(options['blasts'])
        self.build_legend()
        self.build_captions()
        self.build_tracks()
        self.build_cgview()

    def initialize_cgview(self):
        return {
            'version': self.VERSION,
            'created': time.strftime("%Y-%m-%d %H:%M:%S"),
            'id': self.map_id,
            'name': self.map_name,
            'geneticCode': 11,
            'settings': {},
            'backbone': {},
            'ruler': {},
            'dividers': {},
            'annotation': {},
            'sequence': {},
            'captions': [],
            'legend': {},
            'features': [],
            'tracks': []
        }

    def read_sequence(self, path):
        # Determine Sequence file type
        self.seq_type = "genbank"

        # Extract sequence
        sequence_num = 0
        sequence_length = 0
        self.contigs = []
        contig_names = []

        if self.seq_type == "raw":
            with open(path, "r") as f:
                sequence = "".join(line.strip() for line in f if line.strip())
            sequence_num += 1
            sequence_length = len(sequence)
            self.contigs.append(SeqRecord(id="Sequence", seq=sequence))
        else:
            for i, record in enumerate(SeqIO.parse(path, self.seq_type)):
                seq = str(record.seq)
                if not seq:
                    continue
                sequence_num += 1

                contig_temp_name = record.id
                contig_temp_name = contig_temp_name.replace('|', '_')
                contig_name = self.unique_name(contig_temp_name, contig_names)
                contig = {
                    'name': contig_name,
                    'orientation': '+',
                    'length': len(seq),
                    'seq': seq.upper()
                }
                # contig = SeqRecord(id=contig_name, seq=seq.upper())
                sequence_length += len(seq)
                contig_names.append(record)

                self.extract_features(record, contig)

                if i == 0:
                    # self.map_name = record.description
                    if not self.map_name:
                        self.map_name = record.name
                    self.cgview['name'] = self.map_name
                self.contigs.append(contig)
            # print("")

        # print(f"Extracted Sequences: {sequence_num}, Length: {sequence_length} bp, Type: {self.seq_type}")
        # print(f"Extracted Features: {len(self.features)}")

    def detect_filetype(self, path):
        first_line = ""
        with open(path, 'r') as file:
            for line in file:
                if line.strip():
                    first_line = line
                    break

        if re.match(r'^LOCUS\s+', first_line):
            return 'genbank'
        elif re.match(r'^ID\s+', first_line):
            return 'embl'
        elif first_line.startswith('>'):
            return 'fasta'
        else:
            return 'raw'

    def extract_features(self, seq_record, contig=None):
        if self.seq_type not in ('genbank', 'embl'):
            return

        features_to_skip = ['source', 'gene', 'exon']

        for feature in seq_record.features:
            feature_type = feature.type
            if feature_type in features_to_skip or not feature.location:
                continue

            location = feature.location
            start = location.start + 1  # Convert to 1-based index
            end = location.end
            strand = location.strand

            qualifiers = feature.qualifiers
            name = qualifiers.get('gene') or qualifiers.get('locus_tag') or qualifiers.get('note') or qualifiers.get('product') or feature_type
            name = name[0] if name else feature_type  # Handle cases where name is a list
            name = name.encode('utf-8').decode('utf-8')

            codon_start = qualifiers.get('codon_start')
            transl_table = qualifiers.get('transl_table')

            feature_data = {
                'type': feature_type,
                'name': name,
                'start': start,
                'stop': end,
                'strand': strand,
                'source': f'{self.seq_type}-features'
            }

            if contig:
                feature_data['contig'] = contig['name']

            if codon_start and codon_start != '1':
                feature_data['codonStart'] = int(codon_start[0])

            if feature_type == 'CDS':
                genetic_code = int(transl_table[0]) if transl_table else 1
                feature_data['geneticCode'] = genetic_code

            self.features.append(feature_data)

    def build_legend(self):
        config_items = {}
        default_legend_name = None

        if 'legend' in self.config and 'items' in self.config['legend']:
            for item in self.config['legend']['items']:
                config_items[item['name']] = item
            default_legend_name = self.config['legend'].get('default')

        for feature in self.features:
            if feature['type'] in config_items:
                feature['legend'] = config_items[feature['type']]['name']
            elif default_legend_name:
                feature['legend'] = config_items[default_legend_name]['name']
            else:
                feature['legend'] = feature['type']

        feature_legend_names = list(set(config_items.keys()) & set(feature['legend'] for feature in self.features))
        items = [config_items[name] for name in feature_legend_names]

        self.cgview['legend']['items'] = items

    def build_captions(self):
        self.captions = []
        config_captions = self.config.get('captions', [])

        for caption in config_captions:
            if caption['name'].lower() == 'title' and self.map_title() != "":
                caption['name'] = self.map_title()
            self.captions.append(caption)

    def read_blasts(self, paths):
        for i, path in enumerate(paths):
            num = i + 1
            print(f'Creating features for BLAST {num} results...')
            with open(path, 'r') as file:
                reader = csv.reader(file, delimiter='\t')
                for row in reader:
                    query_id = row[0]
                    match_id = row[1]
                    start = int(row[6])
                    stop = int(row[7])
                    strand = 1
                    offset = 0

                    match = re.match(r'([^\t]+)_start=(\d+);end=(\d+)', query_id)
                    if match:
                        offset = int(match.group(2)) - 1

                    meta = {
                        'identity': float(row[2]),
                        'mismatches': int(row[4]),
                        'evalue': float(row[10]),
                        'score': int(row[11])
                    }

                    if start > stop:
                        start, stop = stop, start
                        strand = -1

                    self.features.append({
                        'type': 'blast',
                        'meta': meta,
                        'start': offset + start,
                        'stop': offset + stop,
                        'strand': strand,
                        'score': round(meta['identity'] / 100, 3),
                        'source': f'blast_{num}'
                    })

            self.blast_tracks.append({
                'name': f'blast_{num}',
                'position': 'inside',
                'separateFeaturesBy': 'none',
                'dataType': 'feature',
                'dataMethod': 'source',
                'dataKeys': f'blast_{num}'
            })

    def parse_query_id(self, id):
        query = {'start': 1}
        match = re.match(r'([^\t]+)_start=(\d+);end=(\d+)', id)
        if match:
            query['start'] = int(match.group(2))
            query['end'] = int(match.group(3))
        return query

    def build_tracks(self):
        self.tracks.append({
            'name': 'Features',
            'separateFeaturesBy': 'strand',
            'position': 'both',
            'dataType': 'feature',
            'dataMethod': 'source',
            'dataKeys': f'{self.seq_type}-features'
        })

        if self.blast_tracks:
            self.tracks += self.blast_tracks

        if self.options.get('analysis_path'):
            self.tracks.append({
                'name': 'Analysis',
                'position': 'inside',
                'dataType': 'plot',
                'dataMethod': 'source',
                'dataKeys': 'analysis'
            })

    def build_genetic_code(self):
        codes = defaultdict(int)
        cds_count = 0

        for feature in self.features:
            code = feature.get('geneticCode', self.cgview['geneticCode'])
            codes[code] += 1
            if feature['type'] == 'CDS':
                cds_count += 1

        common_genetic_code = max(codes, key=codes.get)
        print(f"Most Common Genetic Code: {common_genetic_code} (Count: {codes[common_genetic_code]}/{cds_count} CDS)")

        for feature in self.features:
            if feature.get('geneticCode') == common_genetic_code:
                del feature['geneticCode']

    def build_cgview(self):
        print("Creating CGView JSON")
        if self.debug:
            self.cgview['sequence']['seq'] = "SEQUENCE WOULD GO HERE"
            self.cgview['features'] += self.features[1:6]
        else:
            self.cgview['sequence']['contigs'] = self.contigs
            self.cgview['features'] += self.features

        self.cgview['tracks'] = self.tracks + self.cgview['tracks']

        if self.plots:
            self.cgview['plots'] = self.plots

        if self.common_genetic_code:
            self.cgview['geneticCode'] = self.common_genetic_code

        self.cgview['captions'] += self.captions

    def map_title(self):
        if self.map_name:
            return self.map_name
        elif self.seq_type == 'raw':
            return ''
        else:
            return self.seq_record.description

    def to_json(self):
        return json.dumps({'cgview': self.cgview})

    def write_json(self, path):
        with open(path, 'w') as file:
            file.write(self.to_json())

    def unique_name(self, name, all_names):
        if name in all_names:
            count = 2
            new_name = ''
            while True:
                new_name = f"{name}-{count}"
                if new_name not in all_names:
                    break
                count += 1
            return new_name
        return name

    def scale_score(self, score, min_score, max_score):
        return (score - min_score) / (max_score - min_score)


def loadjson(file_path):
    with open(file_path, "r") as json_file:
        data = json.load(json_file)
    return data

@app.route("/")
def index():
    cgview_options = {
        'config': None,
        'map_id': int.from_bytes(os.urandom(4), byteorder="big"),
        'map_name': "coba",
        'contigs': None,
    }
    data = CGViewBuilder(genome_file, options=cgview_options)
    return render_template("index.html", data=data.to_json())

if __name__ == "__main__":
    app.run(debug=True)