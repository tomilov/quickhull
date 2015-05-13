#! /usr/bin/octave -qf

pkg load geometry

lo = -1
hi = 1

facet = unifrnd(lo, hi, 3);
line = [lo lo lo hi - lo hi - lo hi - lo];#unifrnd(lo, hi, 1, 6);#
plane = [facet(1, :) (facet(2, :) - facet(1, :)) (facet(3, :) - facet(1, :))];
point = intersectLinePlane(line, plane);
barycentric = facet'\point';
in_range = all(0 <= barycentric) && all(barycentric <= 1);
coeffs_sum = sum(barycentric);
assert(sumsq(facet' * barycentric - point') < eps);

printf("set label 1 at %f, %f, %f;\n", hi, hi, hi);
printf("set label 1 '%f %f %f';\n", barycentric(1), barycentric(2), barycentric(3));

printf("set xrange [%f:%f];\n", lo, hi);
printf("set yrange [%f:%f];\n", lo, hi);
printf("set zrange [%f:%f];\n", lo, hi);
printf("set view equal xyz;\n");
printf("set xyplane at 0;\n");
if (in_range)
	printf("set title 'hit %f';\n", coeffs_sum);
else
	printf("set title 'miss %f';\n", coeffs_sum);
endif
printf("splot '-' with lines notitle, '-' with lines notitle, '-' with points notitle;\n");

printf("%f ", facet(1, :));
printf("\n");

printf("%f ", facet(2, :));
printf("\n");

printf("%f ", facet(3, :));
printf("\n");

printf("%f ", facet(1, :));
printf("\n");

printf("e\n");


printf("%f ", line(1:3));
printf("\n");
printf("%f ", line(4:6) + line(1:3));
printf("\n");

printf("e\n");


printf("%f ", point(:));
printf("\n");

printf("e\n");

#facet = [facet; facet(1, :)];
#plot3(facet(:, 1), facet(:, 2), facet(:, 3));
