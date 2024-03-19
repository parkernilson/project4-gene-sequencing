#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == "PYQT5":
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == "PYQT4":
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == "PYQT6":
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception("Unsupported Version of PyQt: {}".format(PYQT_VER))

import random
import math

import Proj4GUI

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

BANDED_D = 3


def get_in_bounds(seq1, seq2, banded):
	def in_bounds(row, col):
		if row < 0 or row > len(seq1) or col < 0 or col > len(seq2):
			return False
		if banded == True:
			if col < row - BANDED_D or col > row + BANDED_D:
				return False
			if row < col - BANDED_D or row > row + BANDED_D:
				return False
		return True

	return in_bounds


class RowCol:
	def __init__(self, row, col):
		self.row = row
		self.col = col


class GeneSequencing:

	def __init__(self):
		pass

	def align(self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		seq1_truncated = seq1[:align_length]
		seq2_truncated = seq2[:align_length]

		in_bounds = get_in_bounds(seq1_truncated, seq2_truncated, banded)

		def to_banded(row, col):
			if banded and row > BANDED_D:
				return (row, col - (row - BANDED_D))
			else:
				return (row, col)

		if banded and abs(len(seq1_truncated) - len(seq2_truncated)) > BANDED_D:
			return {
				"align_cost": float("inf"),
				"seqi_first100": "No Alignment Possible.",
				"seqj_first100": "No Alignment Possible.",
			}

		EDIT_DISTANCE = [
			[
				row * 5 if col == 0 else col * 5 if row == 0 else float("inf")
				for col in range(len(seq2_truncated) + 1)
				if in_bounds(row, col)
			]
			for row in range(len(seq1_truncated) + 1)
		]

		PREV = {
			(row, col): (row - 1, col) if col == 0 else (row, col - 1)
			for (row, col) in [(i, 0) for i in range(1, len(seq1_truncated) + 1)]
			+ [(0, j) for j in range(1, len(seq2_truncated) + 1)]
		}

		for row in range(1, len(seq1_truncated) + 1):
			for col in range(1, len(seq2_truncated) + 1):
				if not in_bounds(row, col):
					continue

				(row_i, col_i) = to_banded(row, col) if banded == True else (row, col)

				# prefer left, then top, then diagonal

				if in_bounds(row - 1, col - 1):
					(prev_row_diag, prev_col_diag) = to_banded(row - 1, col - 1) if banded == True else (row - 1, col - 1)
					prev_diag_dist = EDIT_DISTANCE[prev_row_diag][prev_col_diag]
					potential_dist_diag = (
						prev_diag_dist - 3
						if seq1_truncated[row - 1] == seq2_truncated[col - 1]
						else prev_diag_dist + 1
					)
					EDIT_DISTANCE[row_i][col_i] = potential_dist_diag
					PREV[(row, col)] = (row - 1, col - 1)

				(prev_row_top, prev_col_top) = to_banded(row - 1, col) if banded == True else (row - 1, col)
				if in_bounds(row - 1, col):
					prev_top_dist = EDIT_DISTANCE[prev_row_top][prev_col_top]
					potential_dist_top = prev_top_dist + 5
					if potential_dist_top <= EDIT_DISTANCE[row_i][col_i]:
						EDIT_DISTANCE[row_i][col_i] = potential_dist_top
						PREV[(row, col)] = (row - 1, col)

				(prev_row_left, prev_col_left) = to_banded(row, col - 1) if banded == True else (row, col - 1)
				if in_bounds(row, col - 1):
					prev_left_dist = EDIT_DISTANCE[prev_row_left][prev_col_left]
					potential_dist_left = prev_left_dist + 5
					if potential_dist_left <= EDIT_DISTANCE[row_i][col_i]:
						EDIT_DISTANCE[row_i][col_i] = potential_dist_left
						PREV[(row, col)] = (row, col - 1)


		seq_1_aligned_backwards = ""
		seq_2_aligned_backwards = ""
		(cur_row, cur_col) = (len(seq1_truncated), len(seq2_truncated))
		while cur_row > 0 or cur_col > 0:
			(prev_row, prev_col) = PREV[(cur_row, cur_col)]
			# if diag, match or sub
			if prev_row < cur_row and prev_col < cur_col:
				seq_1_aligned_backwards += seq1_truncated[cur_row - 1]
				seq_2_aligned_backwards += seq2_truncated[cur_col - 1]
			# if top, gap for seq 2
			elif prev_row < cur_row:
				seq_1_aligned_backwards += seq1_truncated[cur_row - 1]
				seq_2_aligned_backwards += "-"
			# if left, gap for seq 1
			elif prev_col < cur_col:
				seq_1_aligned_backwards += "-"
				seq_2_aligned_backwards += seq2_truncated[cur_col - 1]
			(cur_row, cur_col) = (prev_row, prev_col)

		seq_1_aligned = seq_1_aligned_backwards[::-1]
		seq_2_aligned = seq_2_aligned_backwards[::-1]

		(end_row, end_col) = (
			to_banded(len(seq1_truncated), len(seq2_truncated))
			if banded
			else (len(seq1_truncated), len(seq2_truncated))
		)
		score = EDIT_DISTANCE[end_row][end_col]
		alignment1 = seq_1_aligned
		alignment2 = seq_2_aligned

		return {
			"align_cost": score,
			"seqi_first100": alignment1[:100],
			"seqj_first100": alignment2[:100],
		}


def loadSequencesFromFile():
	FILENAME = "genomes.txt"
	raw = open(FILENAME, "r").readlines()
	sequences = {}

	i = 0
	cur_id = ""
	cur_str = ""
	for liner in raw:
		line = liner.strip()
		if "#" in line:
			if len(cur_id) > 0:
				sequences[i] = (i, cur_id, cur_str)
				cur_id = ""
				cur_str = ""
				i += 1
			parts = line.split("#")
			cur_id = parts[0]
			cur_str += parts[1]
		else:
			cur_str += line
	if len(cur_str) > 0 or len(cur_id) > 0:
		sequences[i] = (i, cur_id, cur_str)
	return sequences


def testme():
	seqs = loadSequencesFromFile()
	seqs = list(seqs[i][2] for i in seqs)

	seq1 = seqs[8]
	seq2 = seqs[9]

	# seq1_100 = seq1[:10]
	# seq2_100 = seq2[:10]

	gs = GeneSequencing()
	result = gs.align(seq1, seq2, True, 3000)
	print(result)


def test_in_bounds_banded():
	seq1 = "AGCATGC"
	seq2 = "ACAATCC"
	in_bounds = get_in_bounds(seq1, seq2, True)
	assert in_bounds(0, 0) == True
	assert in_bounds(0, 3) == True
	assert in_bounds(0, 4) == False
	assert in_bounds(3, 0) == True
	assert in_bounds(4, 0) == False
	assert in_bounds(3, 7) == False
	assert in_bounds(7, 3) == False
	assert in_bounds(7, 7) == True
	assert in_bounds(8, 7) == False
	assert in_bounds(7, 8) == False


testme()
# test_in_bounds_banded()
