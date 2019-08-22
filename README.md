## proton_beam
Code for "Optimal Targeting of a Tumor through Proton Beam Therapy" with Kiran Pant, currently under review.

Brief description of the files:

1D_beam.py: Generates a series of 1D Bragg-Kleeman curves; Fig. 1 in the paper.

2D_beam.py: Generates a figure from the output of 2D_beam_iteration.py; Fig. 3 in the paper.

2D_beam_iteration.py: Runs simulations for 2D energy deposition of proton beams targeting a simulated tumor; stores data to file.

3D_beam.py: Runs simulation for 3D energy deposition of a proton beam; stores data to file.

3D_beam_plot.py: Generates a figure from the output of 3D_beam.py; Fig. 4 in the paper.

objective_function.py: Compares data generated in 2D_beam_iteration.py; Fig. 2 in the paper.
