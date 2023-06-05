from tempfile import TemporaryDirectory
from werkzeug.utils import secure_filename
import os
from flask import flash
import pandas as pd
from io import BytesIO
from typing import List
import zipfile
import string
from Bio.Seq import Seq
from tracetrack.entities.alignment import Alignment
from tracetrack.entities.aligned import AlignedTrace
from tracetrack.entities.record import TraceSeqRecord
from tracetrack.entities.errors import UnresolvedConflictError
from tracetrack.entities.scheduler import scheduler
from tracetrack.tasks import align_sequences_celery, Settings

ALLOWED_EXTENSIONS = {'zip', 'ab1'}
ALLOWED_REF_EXTENSIONS = {'csv', 'xlsx', 'fasta', 'fa', 'gb', 'gbk', 'genbank'}


def unzip_and_get_sequences(dir: str, mixed_fraction: float) -> List[TraceSeqRecord]:
    """
    Process all files in given directory, return list of TraceSeqRecord objects.
    :param dir: str
    :param threshold: int, score threshold for a position to be taken into account
    :param end_threshold: int, score threshold for low-quality end trimming
    :return: list of TraceSeqRecord objects
    """
    filelist = load_all_files(dir)  # list of file names ending with .ab1
    seqlist = []
    for f in filelist:
        try:
            rec = TraceSeqRecord.read(f, mixed_fraction)
            seqlist.append(rec)
        except ValueError as error:
            flash(str(error))
        except IOError as e:
            flash(f"File {os.path.basename(f)} is in incorrect format.")
    return seqlist


# called by unzip_and_get_sequences
def load_all_files(dir):
    """
    Check all files in directory, unzip zip files, return a list of .ab1 files. Skip invalid files and flash a warning message.
    :param dir: str
    :return: list of file names (str)
    """
    paths = []
    for name in os.listdir(dir):
        if name.lower().endswith('.zip'):
            with zipfile.ZipFile(os.path.join(dir, name), 'r') as zip_ref:
                zip_ref.extractall(dir)
            os.remove(os.path.join(dir, name))
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.lower().endswith('.ab1'):
                paths.append(os.path.join(root, file))
            else:
                flash("File {} has an invalid extension.".format(file))
    # paths.sort()
    return paths


def is_allowed_file(filename, allowed_ext):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in allowed_ext


def save_trace_file(file, tracedir):
    if not is_allowed_file(file.filename, ALLOWED_EXTENSIONS):
        extensions_str = '/'.join(ALLOWED_EXTENSIONS)
        flash(f'Trace file "{file.filename}" type is not recognized, skipping file. '
              f'Please provide trace files in {extensions_str} format.')
        return None
    filename = secure_filename(file.filename)
    file.save(os.path.join(tracedir, filename))
    return filename


def save_trace_files(file_list, tracedir):
    first_list = [save_trace_file(file_name, tracedir) for file_name in file_list]
    return [item for item in first_list if item is not None]


class FileError(Exception):
    def __init__(self, message):
        self.message = message


def get_trace_files(files, num):
    tmp = TemporaryDirectory()

    saved_files = save_trace_files(files, tmp.name)
    if len(saved_files) == 0:
        raise FileError("No selected files for population number {}.".format(num))
    return tmp


def get_ref_input(references):
    if not is_allowed_file(references.filename, ALLOWED_REF_EXTENSIONS):
        extensions_str = '/'.join(ALLOWED_REF_EXTENSIONS)
        raise FileError(f'Reference sheet extension not supported: {references.filename}. '
                        f'Please provide a table in {extensions_str} format.')
    tmp = TemporaryDirectory()
    reffile = secure_filename(references.filename)
    references.save(os.path.join(tmp.name, reffile))
    return tmp, reffile


MUTATION_STATUS_COLORS = {
    'None': '#baffbb',
    'Silent': '#a2cafc',
    'Missense': '#f5a590',
    'Nonsense': '#c19dfa',
    'Frameshift': '#f5d99a'
}


class ResultsError(Exception):
    def __init__(self, message):
        self.message = message


def create_table(alignments: List[Alignment]) -> pd.DataFrame:
    entries = []
    for alignment in alignments:
        ref_id = alignment.aligned_reference.record.id
        filenames = ', '.join(alignment.get_trace_ids())
        mut_type = 'None'
        if alignment.mutations > 0:
            mut_type = 'Silent'
        if alignment.mis_mutations > 0:
            mut_type = 'Missense'
        if alignment.nons_mutations > 0:
            mut_type = 'Nonsense'
        if alignment.frameshifts > 0:
            mut_type = "Frameshift"
        num_reads = len(alignment.aligned_traces)
        entries.append([
            ref_id,
            alignment.population,
            alignment.perc_coverage,
            alignment.perc_identity,
            mut_type,
            alignment.silent_mutations,
            alignment.mis_mutations,
            alignment.nons_mutations,
            alignment.frameshifts,
            num_reads,
            filenames
        ])
    table = pd.DataFrame(entries, columns=['Reference ID', 'Group', '% Coverage', '% Identity', 'Mutation Type',
                                           'Silent Mutations', 'Missense Mutations', 'Nonsense Mutations',
                                           'Frameshift positions',
                                           'Number of reads', 'File names'])

    return table


def number_prefix(string):
    i = 0
    while i < len(string) and string[i].isdigit():
        i += 1
    return i


def shorten_sheet_names(names_orig, max_length=27, max_iter=25):
    suffixes = [0] * len(names_orig)
    names_final = [name[:max_length] for name in names_orig]
    iteration = 0
    while len(set(names_final)) < len(names_final) and iteration < max_iter:
        iteration += 1
        used = {}
        for i, name in enumerate(names_final):
            if names_final.count(name) > 1:
                number = used.get(name, 0) + 1
                used[name] = number
                suffixes[i] = number
        for i in range(len(names_orig)):
            suffix = f" ({str(suffixes[i])})" if suffixes[i] > 0 else ""
            if max_length - len(suffix) < 0:
                raise ValueError(f"Suffix of sheet name {suffix} is too long!")
            if len(names_orig[i]) > max_length - len(suffix):
                suffix = "..." + suffix
            names_final[i] = f"{names_orig[i][:max_length - len(suffix)]}{suffix}"
    if len(set(names_final)) < len(names_final):
        collisions = set([name for name in names_final if names_final.count(name) > 1])
        raise UnresolvedConflictError(f"The colliding sheet names could not be resolved. Collisions: {collisions}")
    return names_final


def create_spreadsheet(alignments: List[Alignment], results_table=None, index=None):
    if results_table is None:
        results_table = create_table(alignments)
    output = BytesIO()
    writer = pd.ExcelWriter(output)

    if results_table is None:
        raise ResultsError("Nothing to download.")

    # write result table to spreadsheet
    results_table.to_excel(writer, startrow=0, merge_cells=False, sheet_name="Results")
    workbook = writer.book
    wrap_format = workbook.add_format({'text_wrap': True})

    # Formatting settings for the spreadsheet
    # background colors for "Mutation" column
    none_format = workbook.add_format({'bg_color': MUTATION_STATUS_COLORS['None']})
    mis_format = workbook.add_format({'bg_color': MUTATION_STATUS_COLORS['Missense']})
    silent_format = workbook.add_format({'bg_color': MUTATION_STATUS_COLORS['Silent']})
    nons_format = workbook.add_format({'bg_color': MUTATION_STATUS_COLORS['Nonsense']})
    frameshift_format = workbook.add_format({'bg_color': MUTATION_STATUS_COLORS['Frameshift']})

    # formats for alignment positions
    no_mut = workbook.add_format({'font_color': 'black'})
    no_mut_dark = workbook.add_format({'font_color': 'black', 'bg_color': '#e8e8e8'})
    silent_mut = workbook.add_format({'bold': True, 'font_color': '#0070fa', 'underline': True})
    silent_mut_dark = workbook.add_format(
        {'bold': True, 'font_color': '#0070fa', 'underline': True, 'bg_color': '#e8e8e8'})
    mis_mut = workbook.add_format({'bold': True, 'font_color': '#fc5e42', 'underline': True})
    mis_mut_dark = workbook.add_format(
        {'bold': True, 'font_color': '#fc5e42', 'underline': True, 'bg_color': '#e8e8e8'})
    nons_mut = workbook.add_format({'bold': True, 'font_color': '#6d23ad', 'underline': True, 'italic': True})
    nons_mut_dark = workbook.add_format(
        {'bold': True, 'font_color': '#6d23ad', 'underline': True, 'italic': True, 'bg_color': '#e8e8e8'})
    frameshift = workbook.add_format({'bold': True, 'font_color': '#fcbd28', 'underline': True, 'italic': True})
    frameshift_dark = workbook.add_format(
        {'bold': True, 'font_color': '#fcbd28', 'underline': True, 'italic': True, 'bg_color': '#e8e8e8'})
    formats = {'no_mut': no_mut, 'mis_mut': mis_mut, 'silent_mut': silent_mut, 'nons_mut': nons_mut,
               'frameshift': frameshift}
    formats_dark = {'no_mut': no_mut_dark, 'mis_mut': mis_mut_dark, 'silent_mut': silent_mut_dark,
                    'nons_mut': nons_mut_dark, 'frameshift': frameshift_dark}

    rows = len(results_table) + 1

    worksheet = writer.sheets["Results"]

    # Conditional formatting based on coverage and identity
    conditional_format_by_bins(workbook, worksheet, rows)

    # Conditional formatting based on mutation type
    worksheet.conditional_format('F2:F{}'.format(rows), {'type': 'text',
                                                         'criteria': 'containing',
                                                         'value': 'Nonsense',
                                                         'format': nons_format})
    worksheet.conditional_format('F2:F{}'.format(rows), {'type': 'text',
                                                         'criteria': 'containing',
                                                         'value': 'Missense',
                                                         'format': mis_format})
    worksheet.conditional_format('F2:F{}'.format(rows), {'type': 'text',
                                                         'criteria': 'containing',
                                                         'value': 'Silent',
                                                         'format': silent_format})
    worksheet.conditional_format('F2:F{}'.format(rows), {'type': 'text',
                                                         'criteria': 'containing',
                                                         'value': 'None',
                                                         'format': none_format})
    worksheet.conditional_format('F2:F{}'.format(rows), {'type': 'text',
                                                         'criteria': 'containing',
                                                         'value': 'Frameshift',
                                                         'format': frameshift_format})

    worksheet.set_column('A:K', 15)

    # formats for alignment display
    rotated = workbook.add_format({'bottom': 2})
    rotated.set_rotation(90)
    bottom = workbook.add_format({'bottom': 2})
    bold = workbook.add_format({'bold': True})
    bold_dark = workbook.add_format({'bold': True, 'bg_color': '#e8e8e8'})
    bold_dark_bottom = workbook.add_format({'bold': True, 'bg_color': '#e8e8e8', 'bottom': 2})
    bold_bottom = workbook.add_format({'bold': True, 'bottom': 2})
    dark = workbook.add_format({'bg_color': '#e8e8e8'})
    dark_bottom = workbook.add_format({'bg_color': '#e8e8e8', 'bottom': 2})

    worksheet_res = worksheet

    sheet_names = shorten_sheet_names([al.ref_id for al in alignments])
    # Iterate through all alignments, create sheet for each
    for num, alignment in enumerate(alignments):
        name = sheet_names[num]
        worksheet = workbook.add_worksheet(name + " Seq")
        worksheet_mut = workbook.add_worksheet(name + " Mut")
        # create link from results table to alignment sheet
        if index is not None:
            worksheet_res.write_url(f'B{index + 2}', f'internal:\'{name} Seq\'!A1', string=alignment.ref_id)
        else:
            worksheet_res.write_url(f'B{num + 2}', f'internal:\'{name} Seq\'!A1', string=alignment.ref_id)

        # write alignment sequences: reference and mutated (result), DNA and AA
        row = 0
        col = 0
        worksheet.set_column(1, len(alignment.positions) + 10, 1.8)
        worksheet.set_column(0, 0, 14)
        worksheet.set_row(0, None, rotated)
        worksheet.write(1, col, 'Reference DNA', bold)
        worksheet.write(2, col, 'Reference AA', bold_bottom)
        worksheet.write(3, col, 'Consensus DNA')
        worksheet.write(4, col, 'Consensus AA', bottom)

        col += 6
        # is Kozak sequence present?
        if alignment.kozak_count:
            worksheet.write(3, 1, 'Kozak present')
        else:
            worksheet.write(3, 1, 'Kozak missing', mis_mut)

        col += 1
        col_beg = col
        codon_count_ref = 3
        codon_ref = []
        codon_count_cons = 3
        codon_cons = []
        dark_ref = True
        dark_cons = True
        for i, pos in enumerate(alignment.positions):
            style_ref = (bold, bold_dark)[dark_ref]
            style_ref_aa = bold_dark_bottom if dark_ref else bold_bottom
            style_cons_aa = dark_bottom if dark_cons else bottom

            # position numbering
            worksheet.write(row, col, alignment.position_to_display(i), rotated)
            # DNA ref
            worksheet.write(row + 1, col, pos.ref, style_ref)
            # AA ref
            worksheet.write(row + 2, col, "", style_ref_aa)
            if pos.ref != '-' and pos.ref != '?':
                codon_count_ref = (codon_count_ref + 1) % 3
                codon_ref.append(pos.ref)
            if codon_count_ref == 0:
                aa = str(Seq("".join(codon_ref)).translate())
                codon_ref = []
                worksheet.write(row + 2, col - 1, aa, style_ref_aa)
                dark_ref = not dark_ref
            # DNA consensus
            if dark_cons:
                worksheet.write(row + 3, col, pos.result, formats_dark[alignment.mutation_type(i)])
            else:
                worksheet.write(row + 3, col, pos.result, formats[alignment.mutation_type(i)])
            # AA consensus
            worksheet.write(row + 4, col, "", style_cons_aa)
            if pos.result != '-' and pos.result != '?':
                codon_count_cons = (codon_count_cons + 1) % 3
                codon_cons.append(pos.result)
            if codon_count_cons == 0:
                aa = str(Seq("".join(codon_cons)).translate())
                codon_cons = []
                worksheet.write(row + 4, col - 1, aa, style_cons_aa)
                dark_cons = not dark_cons
            col += 1

        worksheet.set_row(row + 2, None, bottom)
        worksheet.set_row(row + 4, None, bottom)
        row += 5

        # for each aligned trace file, write sequence under ref and result sequences
        for aligned_trace in alignment.aligned_traces:
            codon_count = 3
            codon = []
            dark_pos = True
            col = col_beg
            worksheet.write(row, 0, aligned_trace.record.id)
            # write leading/Kozak sequence, if present
            if len(aligned_trace.leading_seq) >= 6:
                kozak = aligned_trace.leading_seq[-6:]
                for i in range(6):
                    worksheet.write(row, i + 1, kozak[i])
            for pos in aligned_trace.aligned_seq:
                style = dark if dark_pos else None
                style_aa = dark_bottom if dark_pos else bottom
                worksheet.write(row, col, pos, style)
                worksheet.write(row + 1, col, "", style_aa)
                if pos != '-':
                    codon_count = (codon_count + 1) % 3
                    codon.append(pos)
                if codon_count == 0:
                    aa = str(Seq("".join(codon)).translate())
                    codon = []
                    worksheet.write(row + 1, col - 1, aa, style_aa)
                    dark_pos = not dark_pos
                col += 1
            worksheet.set_row(row + 1, None, bottom)
            # stop codon, if present
            if len(aligned_trace.trailing_seq) >= 3:
                for i in range(3):
                    worksheet.write(row, col, aligned_trace.trailing_seq[i])
            row += 2

        # is a stop codon present?
        style = no_mut if alignment.stop_count else mis_mut
        worksheet.write(1, col, '*', style)

        row = 0
        row_mark = row
        # write positions of all mismatches
        worksheet_mut.write(row, 0, "Mismatches", bold)
        row += 1
        worksheet_mut.write(row, 0, "Position", rotated)
        worksheet_mut.write(row, 1, "Reference", rotated)
        worksheet_mut.write(row, 2, "Result", rotated)
        worksheet_mut.write(row, 3, "Mutation Type", rotated)
        worksheet_mut.write(row, 4, "Amino Acid Change", rotated)
        worksheet_mut.write(row, 5, "Comment", rotated)
        worksheet_mut.set_column(0, 0, 3.8)
        worksheet_mut.set_column(1, 2, 2)
        worksheet_mut.set_column(3, 3, 8.5)
        worksheet_mut.set_column(4, 4, 5)
        for mismatch in alignment.mismatches:
            row += 1
            worksheet_mut.write_url(row, 0, f'internal:\'{name} Seq\'!{get_cell_column(mismatch + 8)}1',
                                    string=str(mismatch + 1))
            worksheet_mut.write(row, 1, alignment.positions[mismatch].ref)
            worksheet_mut.write(row, 2, alignment.positions[mismatch].result)
            if alignment.positions[mismatch].result == "?":
                worksheet_mut.write(row, 5, "(\"?\" is an unknown insertion)")
            worksheet_mut.write(row, 3, alignment.mutation_type(mismatch))
            original_aa = alignment.translate_reference(mismatch)
            mut_aa = alignment.translate_mutated(mismatch)
            worksheet_mut.write(row, 4, f"{original_aa} -> {mut_aa}")
        col = 11
        row = row_mark
        no_coverage = alignment.get_zero_coverage_runs()

        # write locations of all runs with zero coverage
        worksheet_mut.write(row, col, "Zero coverage:", bold)
        row += 1
        for pos in no_coverage:
            worksheet_mut.write(row, col, str(pos))
            row += 1

    writer.close()

    # go back to the beginning of the stream
    output.seek(0)
    return output


color_bins = [
    '#f04646',
    '#ff5353',
    '#ff8039',
    '#ffac33',
    '#fac151',
    '#96c66a',
    '#5fc85f'
]


def conditional_format_by_bins(workbook, worksheet, row_num):
    # background colors for "Coverage" column
    cov_90 = workbook.add_format({'bg_color': color_bins[6]})
    cov_85 = workbook.add_format({'bg_color': color_bins[5]})
    cov_80 = workbook.add_format({'bg_color': color_bins[4]})
    cov_75 = workbook.add_format({'bg_color': color_bins[3]})
    cov_70 = workbook.add_format({'bg_color': color_bins[2]})
    cov_65 = workbook.add_format({'bg_color': color_bins[1]})
    cov_00 = workbook.add_format({'bg_color': color_bins[0]})

    worksheet.conditional_format('D2:D{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 0,
                                                            'maximum': 65,
                                                            'format': cov_00})
    worksheet.conditional_format('D2:D{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 65,
                                                            'maximum': 70,
                                                            'format': cov_65})
    worksheet.conditional_format('D2:D{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 70,
                                                            'maximum': 75,
                                                            'format': cov_70})
    worksheet.conditional_format('D2:D{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 75,
                                                            'maximum': 80,
                                                            'format': cov_75})
    worksheet.conditional_format('D2:D{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 80,
                                                            'maximum': 85,
                                                            'format': cov_80})
    worksheet.conditional_format('D2:D{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 85,
                                                            'maximum': 90,
                                                            'format': cov_85})
    worksheet.conditional_format('D2:D{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 90,
                                                            'maximum': 101,
                                                            'format': cov_90})
    # Identity
    worksheet.conditional_format('E2:E{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 0,
                                                            'maximum': 95,
                                                            'format': cov_00})
    worksheet.conditional_format('E2:E{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 95,
                                                            'maximum': 96,
                                                            'format': cov_65})
    worksheet.conditional_format('E2:E{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 96,
                                                            'maximum': 97,
                                                            'format': cov_70})
    worksheet.conditional_format('E2:E{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 97,
                                                            'maximum': 98,
                                                            'format': cov_75})
    worksheet.conditional_format('E2:E{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 98,
                                                            'maximum': 99,
                                                            'format': cov_80})
    worksheet.conditional_format('E2:E{}'.format(row_num), {'type': 'cell',
                                                            'criteria': 'between',
                                                            'minimum': 99,
                                                            'maximum': 100,
                                                            'format': cov_85})
    worksheet.conditional_format('E2:E{}'.format(row_num), {'type': 'cell',
                                                            'criteria': '>=',
                                                            'value': 100,
                                                            'format': cov_90})


def get_cell_column(index):
    """Convert 1-based index to alphabetical label of excel column"""
    res = ""
    letters = [""] + list(string.ascii_uppercase)
    while index > 0:
        modulo = index % 26
        if modulo == 0:
            modulo = 26
            index -= 1
        res = letters[modulo] + res
        index = index // 26
    return res


def is_best_for_reference(results, row, field):
    reference = row['Reference ID']
    return results[results['Reference ID'] == reference][field].max() == row[field]


def get_color_bin(value, thresholds):
    bin = 0
    for i, threshold in enumerate(thresholds):
        if value >= threshold:
            bin = i + 1
    return bin


def trace_record_dict(trace: AlignedTrace) -> dict:
    return {
        'sequence': str(trace.record.seq),
        'traces': [{'base': base, 'values': trace.record.traces[base]} for base in 'GATC'],
        'locations': trace.record.base_locations,
        'alignedPositions': trace.get_aligned_positions(),
        'id': trace.record.id
    }


def position_dict(alignment: Alignment, position: int) -> dict:
    return {
        'mut': alignment.mutation_type(position),
        'ref': alignment.positions[position].pos_ref,
        'cov': alignment.positions[position].coverage,
        'pos': alignment.position_to_display(position)
    }


def has_extension(filename, exts):
    if filename.split(".")[-1].lower() in exts:
        return True
    return False


def hash_trace_name(pop_number, seq_number):
    return f"{pop_number}_{seq_number}"


def schedule_tasks(sequences, population_names, db, separate, threshold, end_threshold):
    settings = Settings(threshold, end_threshold, separate)
    inputs = []
    inputs_sorted = []
    for population_seqs, population_name in zip(sequences, population_names):
        new_seqlist = []
        warnings = []
        for record in population_seqs:
            new_record = record.filter_sequence_by_quality(threshold, end_threshold)
            if new_record.has_base_above_threshold():
                new_seqlist.append(new_record)
            else:
                warnings.append(
                    f"Skipping trace file {new_record.id} because no bases passed given quality thresholds. Try resubmitting the task with lower thresholds.")

        if separate:
            # each sequence as separate list with one item
            seqlists = [[s] for s in new_seqlist]
        else:
            # single list with all sequences
            seqlists = [new_seqlist]

        for seqlist in seqlists:
            seqs_by_refid = {}
            for record in seqlist:
                if record.reverse:
                    record = record.reverse_complement()
                seqs_by_refid[record.reference] = seqs_by_refid.get(record.reference, []) + [record]
            for ref_id, seqlist in seqs_by_refid.items():
                inputs.append({
                    "seqlist": seqlist, "ref_id": ref_id, "population_name": population_name, "db": db,
                    "warnings_orig": warnings, "settings": settings
                })
                inputs_sorted = sorted(inputs, key=lambda a: (
                    number_prefix(a["ref_id"]) > 0, int(a["ref_id"][:number_prefix(a["ref_id"])])
                    if number_prefix(a["ref_id"]) > 0 else a["ref_id"]
                ))
                warnings = []
    task_id = scheduler.schedule_tasks(
        align_sequences_celery,
        inputs=inputs_sorted
    )
    return task_id
