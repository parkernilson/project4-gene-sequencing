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

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


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

		EDIT_DISTANCE = [
			[
				row * 5 if col == 0 else col * 5 if row == 0 else math.inf
				for col in range(len(seq2) + 1)
			]
			for row in range(len(seq1) + 1)
		]

		PREV = {
			(row, col): (row - 1, col) if col == 0 else (row, col - 1)
			for (row, col) in [(i, 0) for i in range(1, len(seq1) + 1)]
			+ [(0, j) for j in range(1, len(seq2) + 1)]
		}

		for row in range(1, len(seq1) + 1):
			for col in range(1, len(seq2) + 1):
				# prefer left, then top, then diagonal

				prev_diag_dist = EDIT_DISTANCE[row - 1][col - 1]
				potential_dist_diag = (
					prev_diag_dist - 3
					if seq1[row - 1] == seq2[col - 1]
					else prev_diag_dist + 1
				)
				EDIT_DISTANCE[row][col] = potential_dist_diag
				PREV[(row, col)] = (row - 1, col - 1)

				prev_top_dist = EDIT_DISTANCE[row - 1][col]
				potential_dist_top = prev_top_dist + 5
				if potential_dist_top <= EDIT_DISTANCE[row][col]:
					EDIT_DISTANCE[row][col] = potential_dist_top
					PREV[(row, col)] = (row - 1, col)

				prev_left_dist = EDIT_DISTANCE[row][col - 1]
				potential_dist_left = prev_left_dist + 5
				if potential_dist_left <= EDIT_DISTANCE[row][col]:
					EDIT_DISTANCE[row][col] = potential_dist_left
					PREV[(row, col)] = (row, col - 1)

		seq_1_aligned_backwards = ""
		seq_2_aligned_backwards = ""
		(cur_row, cur_col) = (len(seq1), len(seq2))
		while cur_row > 0 or cur_col > 0:
			(prev_row, prev_col) = PREV[(cur_row, cur_col)]
			# if diag, match or sub
			if prev_row < cur_row and prev_col < cur_col:
				seq_1_aligned_backwards += seq1[cur_row - 1]
				seq_2_aligned_backwards += seq2[cur_col - 1]
			# if top, gap for seq 2
			elif prev_row < cur_row:
				seq_1_aligned_backwards += seq1[cur_row - 1]
				seq_2_aligned_backwards += "-"
			# if left, gap for seq 1
			elif prev_col < cur_col:
				seq_1_aligned_backwards += "-"
				seq_2_aligned_backwards += seq2[cur_col - 1]
			(cur_row, cur_col) = (prev_row, prev_col)

		seq_1_aligned = seq_1_aligned_backwards[::-1]
		seq_2_aligned = seq_2_aligned_backwards[::-1]

		score = EDIT_DISTANCE[len(seq1)][len(seq2)]
		alignment1 = f"""{seq_1_aligned} 
			DEBUG:({len(seq1)} 
			chars,align_len={align_length}{',BANDED' if banded else ''})
			"""
		alignment2 = f"""{seq_2_aligned} 
			DEBUG:({len(seq2)} 
			chars,align_len={align_length}{',BANDED' if banded else ''})
			"""

		return {
			"align_cost": score,
			"seqi_first100": alignment1,
			"seqj_first100": alignment2,
		}


def test_me():
	seq1 = "polynomial"
	seq2 = "exponential"
	gs = GeneSequencing()
	return gs.align(seq1, seq2, False, None)


result = test_me()
print(result)