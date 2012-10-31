# matplotlib-based plotting script for Lab 5 data
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import numpy as np


# scaling results (fixed levels)
nprocs = (1, 8, 27, 64, 125, 216, 343);
np_scaling_l1 = (0.28, 1.05, 1.65, 1.90, 4.57, 2.27, 2.74);
np_scaling_l2 = (6.49, 38.16, 55.29, 64.29, 72.53, 84.18, 94.30);
np_scaling_l3 = (19.33, 115.92, 135.18, 184.58, 192.71, 238.54, 246.95);
np_scaling_l4 = (43.32, 266.44, 296.94, 373.02, 412.97, 490.11, 552.68);
np_scaling_l5 = (82.66, 492.70, 563.73, 663.54, 737.91, 905.84, 911.34);
np_scaling_l6 = (139.68, 845.25, 1177.91, 1082.60, 1254.13, 1462.94, 1521.96);
np_scaling_l7 = (248.57, 1399.96, 1716.65, 1906.60, 1963.85, 2340.41, 2515.80);
p_scaling_l1 = (0.27, 0.85, 1.41, 1.19, 2.37, 2.63, 2.92);
p_scaling_l2 = (10.82, 35.72, 44.15, 43.99, 64.19, 76.27, 81.12);
p_scaling_l3 = (24.46, 83.54, 91.58, 91.46, 130.09, 152.20, 174.85);
p_scaling_l4 = (42.95, 147.29, 163.24, 161.10, 217.19, 269.89, 301.19);
p_scaling_l5 = (69.74, 232.09, 267.64, 262.78, 357.96, 147.19, 457.53);
p_scaling_l6 = (93.96, 344.71, 386.19, 380.37, 529.23, 581.60, 647.26);
p_scaling_l7 = (131.70, 495.90, 553.14, 538.61, 712.36, 821.09, 908.49);
p_scaling_l8 = (290.76, 1067.70, 1282.83, 1338.18, 1568.88, 1774.67, 1982.97);
p_scaling_l9 = (852.99, 2825.14, 4069.76, 4420.10, 5309.89, 5764.71, 7302.69);

levels1 = ('1 level','2 levels','3 levels','4 levels','5 levels','6 levels','7 levels');
levels2 = ('1 level','2 levels','3 levels','4 levels','5 levels','6 levels','7 levels','8 levels','9 levels');


# generate the scaling plots
figure()
loglog(nprocs, np_scaling_l1, 'm-')
loglog(nprocs, np_scaling_l2, 'r-')
loglog(nprocs, np_scaling_l3, 'b-')
loglog(nprocs, np_scaling_l4, 'k-')
loglog(nprocs, np_scaling_l5, 'c-')
loglog(nprocs, np_scaling_l6, 'm--')
loglog(nprocs, np_scaling_l7, 'r--')
xlabel('Processes')
ylabel('Solve Time')
title('Solution Time Weak Scaling (no prec.)')
legend(levels1, loc='lower left', shadow=True)
axis((1, 350, 0.1, 10000))
grid()
savefig('scaling_np.pdf')

figure()
loglog(nprocs, p_scaling_l1, 'm-')
loglog(nprocs, p_scaling_l2, 'r-')
loglog(nprocs, p_scaling_l3, 'b-')
loglog(nprocs, p_scaling_l4, 'k-')
loglog(nprocs, p_scaling_l5, 'c-')
loglog(nprocs, p_scaling_l6, 'm--')
loglog(nprocs, p_scaling_l7, 'r--')
loglog(nprocs, p_scaling_l8, 'b--')
loglog(nprocs, p_scaling_l9, 'k--')
xlabel('Processes')
ylabel('Solve Time')
title('Solution Time Weak Scaling (HG prec.)')
legend(levels2, loc='lower left', shadow=True)
axis((1, 350, -1, 10000))
grid()
savefig('scaling_p.pdf')


# runtime results for fixed processor numbers
grids_p1 = np.array((32768, 32768 + 46656, 32768 + 46656 + 46656, 32768 + 46656 + 46656 + 46656, 32768 + 46656 + 46656 + 46656 + 46656, 32768 + 46656 + 46656 + 46656 + 46656 + 46656, 32768 + 46656 + 46656 + 46656 + 46656 + 46656 + 46656, 32768 + 46656 + 46656 + 46656 + 46656 + 46656 + 46656 + 46656, 32768 + 46656 + 46656 + 46656 + 46656 + 46656 + 46656 + 46656 + 46656), dtype=float);
grids_p8 = np.array((262144, 262144 + 373248, 262144 + 373248 + 373248, 262144 + 373248 + 373248 + 373248, 262144 + 373248 + 373248 + 373248 + 373248, 262144 + 373248 + 373248 + 373248 + 373248 + 373248, 262144 + 373248 + 373248 + 373248 + 373248 + 373248 + 373248, 262144 + 373248 + 373248 + 373248 + 373248 + 373248 + 373248 + 373248, 262144 + 373248 + 373248 + 373248 + 373248 + 373248 + 373248 + 373248 + 373248), dtype=float);
grids_p27 = np.array((884736, 884736 + 1259712, 884736 + 1259712 + 1259712, 884736 + 1259712 + 1259712 + 1259712, 884736 + 1259712 + 1259712 + 1259712 + 1259712, 884736 + 1259712 + 1259712 + 1259712 + 1259712 + 1259712, 884736 + 1259712 + 1259712 + 1259712 + 1259712 + 1259712 + 1259712, 884736 + 1259712 + 1259712 + 1259712 + 1259712 + 1259712 + 1259712 + 1259712, 884736 + 1259712 + 1259712 + 1259712 + 1259712 + 1259712 + 1259712 + 1259712 + 1259712), dtype=float);
grids_p64 = np.array((2097152, 2097152 + 2985984, 2097152 + 2985984 + 2985984, 2097152 + 2985984 + 2985984 + 2985984, 2097152 + 2985984 + 2985984 + 2985984 + 2985984, 2097152 + 2985984 + 2985984 + 2985984 + 2985984 + 2985984, 2097152 + 2985984 + 2985984 + 2985984 + 2985984 + 2985984 + 2985984, 2097152 + 2985984 + 2985984 + 2985984 + 2985984 + 2985984 + 2985984 + 2985984, 2097152 + 2985984 + 2985984 + 2985984 + 2985984 + 2985984 + 2985984 + 2985984 + 2985984), dtype=float);
grids_p125 = np.array((4096000, 4096000 + 5832000, 4096000 + 5832000 + 5832000, 4096000 + 5832000 + 5832000 + 5832000, 4096000 + 5832000 + 5832000 + 5832000 + 5832000, 4096000 + 5832000 + 5832000 + 5832000 + 5832000 + 5832000, 4096000 + 5832000 + 5832000 + 5832000 + 5832000 + 5832000 + 5832000, 4096000 + 5832000 + 5832000 + 5832000 + 5832000 + 5832000 + 5832000 + 5832000, 4096000 + 5832000 + 5832000 + 5832000 + 5832000 + 5832000 + 5832000 + 5832000 + 5832000), dtype=float);
grids_p216 = np.array((7077888, 7077888 + 10077696, 7077888 + 10077696 + 10077696, 7077888 + 10077696 + 10077696 + 10077696, 7077888 + 10077696 + 10077696 + 10077696 + 10077696, 7077888 + 10077696 + 10077696 + 10077696 + 10077696 + 10077696, 7077888 + 10077696 + 10077696 + 10077696 + 10077696 + 10077696 + 10077696, 7077888 + 10077696 + 10077696 + 10077696 + 10077696 + 10077696 + 10077696 + 10077696, 7077888 + 10077696 + 10077696 + 10077696 + 10077696 + 10077696 + 10077696 + 10077696 + 10077696), dtype=float);
grids_p343 = np.array((11239424, 11239424 + 16003008, 11239424 + 16003008 + 16003008, 11239424 + 16003008 + 16003008 + 16003008, 11239424 + 16003008 + 16003008 + 16003008 + 16003008, 11239424 + 16003008 + 16003008 + 16003008 + 16003008 + 16003008, 11239424 + 16003008 + 16003008 + 16003008 + 16003008 + 16003008 + 16003008, 11239424 + 16003008 + 16003008 + 16003008 + 16003008 + 16003008 + 16003008 + 16003008, 11239424 + 16003008 + 16003008 + 16003008 + 16003008 + 16003008 + 16003008 + 16003008 + 16003008), dtype=float);
nlevels1 = (1, 2, 3, 4, 5, 6, 7);
np_runtime_p1   = np.array((0.28, 6.49, 19.33, 43.32, 82.66, 139.68, 248.57)) / grids_p1[0:7];
np_runtime_p8   = np.array((1.05, 38.16, 115.92, 266.44, 492.70, 845.25, 1399.96)) / grids_p8[0:7];
np_runtime_p27  = np.array((1.65, 55.29, 135.18, 296.94, 563.73, 1177.91, 1716.65)) / grids_p27[0:7];
np_runtime_p64  = np.array((1.90, 64.29, 184.58, 373.02, 663.54, 1082.60, 1906.60)) / grids_p64[0:7];
np_runtime_p125 = np.array((4.57, 72.53, 192.71, 412.97, 737.91, 1254.13, 1963.85)) / grids_p125[0:7];
np_runtime_p216 = np.array((2.27, 84.18, 238.54, 490.11, 905.84, 1462.94, 2340.41)) / grids_p216[0:7];
np_runtime_p343 = np.array((2.74, 94.30, 246.95, 552.68, 911.34, 1521.96, 2515.80)) / grids_p343[0:7];

nlevels2 = (1, 2, 3, 4, 5, 6, 7, 8, 9);
p_runtime_p1   = np.array((0.27, 10.82, 24.46, 42.95, 69.74, 93.96, 131.70, 290.76, 852.99)) / grids_p1;
p_runtime_p8   = np.array((0.85, 35.72, 83.54, 147.29, 232.09, 344.71, 495.90, 1067.70, 2825.14)) / grids_p8;
p_runtime_p27  = np.array((1.41, 44.15, 91.58, 163.24, 267.64, 386.19, 553.14, 1282.83, 4069.76)) / grids_p27;
p_runtime_p64  = np.array((1.19, 43.99, 91.46, 161.10, 262.78, 380.37, 538.61, 1338.18, 4420.10)) / grids_p64;
p_runtime_p125 = np.array((2.37, 64.19, 130.09, 217.19, 357.96, 529.23, 712.36, 1568.88, 5309.89)) / grids_p125;
p_runtime_p216 = np.array((2.63, 76.27, 152.20, 269.89, 147.19, 581.60, 821.09, 1774.67, 5764.71)) / grids_p216;
p_runtime_p343 = np.array((2.92, 81.12, 174.85, 301.19, 457.53, 647.26, 908.49, 1982.97, 7302.69)) / grids_p343;

procs = ('1 proc','8 procs','27 procs','64 procs','125 procs','216 procs','343 procs');

# generate the timing plots
figure()
semilogy(nlevels1, np_runtime_p1,   'm-')
semilogy(nlevels1, np_runtime_p8,   'r-')
semilogy(nlevels1, np_runtime_p27,  'b-')
semilogy(nlevels1, np_runtime_p64,  'k-')
semilogy(nlevels1, np_runtime_p125, 'c-')
semilogy(nlevels1, np_runtime_p216, 'm--')
semilogy(nlevels1, np_runtime_p343, 'r--')
xlabel('AMR Levels')
ylabel('Solve Time')
title('Average Solution Time Per Cell With AMR (no prec.)')
legend(procs, loc='lower right', shadow=True)
axis((1, 9, 1e-7, 1e-2))
grid()
savefig('timing_np.pdf')

figure()
semilogy(nlevels2, p_runtime_p1,   'm-')
semilogy(nlevels2, p_runtime_p8,   'r-')
semilogy(nlevels2, p_runtime_p27,  'b-')
semilogy(nlevels2, p_runtime_p64,  'k-')
semilogy(nlevels2, p_runtime_p125, 'c-')
semilogy(nlevels2, p_runtime_p216, 'm--')
semilogy(nlevels2, p_runtime_p343, 'r--')
xlabel('AMR Levels')
ylabel('Solve Time')
title('Average Solution Time Per Cell With AMR (HG prec.)')
legend(procs, loc='upper left', shadow=True)
axis((1, 9, 1e-7, 1e-2))
grid()
savefig('timing_p.pdf')
