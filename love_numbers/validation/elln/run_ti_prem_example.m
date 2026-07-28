work_directory = getenv('ELLN_WORK_DIR');
output_file = getenv('ELLN_OUTPUT');
if isempty(work_directory) || isempty(output_file)
    error('ELLN_WORK_DIR and ELLN_OUTPUT must be set');
end

cd(work_directory);
degrees = [1, 2, 10, 50];
file = fopen(output_file, 'w');
if file == -1
    error('Could not open ELLN output file');
end

for degree = degrees
    [n, h, horizontal, k, mass, gravity] = ELLN( ...
        'EarthMantleTI56.txt', 'EarthCore26.txt', ...
        3, 3, degree, degree, 1);
    fprintf(file, '%d %.17g %.17g %.17g %.17g %.17g\n', ...
        n(1), h(1), horizontal(1), k(1), mass, gravity);
end
fclose(file);
