fileOutd = '/home/mstobbe/FBPIC_Sims/drive_test.txt';
fileOutw = '/home/mstobbe/FBPIC_Sims/witness_test.txt';

conda_environment = 'masontest'; % Replace with your Conda environment name
python_script = '/home/mstobbe/FBPIC_Sims/pwaSim.py'; % Replace with the path to your Python script
arguments = [fileOutd,' ',fileOutw] % Replace with your desired command line arguments

init = 'conda init --all'
activate = 'conda activate masontest'
command = sprintf('conda run --cwd /home/mstobbe/FBPIC_Sims -n masontest python3.11 %s %s', python_script, arguments)
command2 = sprintf('/home/mstobbe/anaconda3/envs/masontest/bin/python %s %s', python_script, arguments)
pip = 'pip install cupy-cuda110'

[status, result] = system(command);
disp(result)

%% venv

fileOutd = '/home/mstobbe/plasmaSims/drive_test.txt';
fileOutw = '/home/mstobbe/plasmaSims/witness_test.txt';


virtualenv_path = '/home/mstobbe/plasmaSims/'; % Replace with the path to your virtual environment
python_script = '/home/mstobbe/FBPIC_Sims/pwaSim.py'; % Replace with the path to your Python script
arguments = [fileOutd,' ',fileOutw] % Replace with your desired command line arguments

activate_script = fullfile(virtualenv_path, 'bin', 'activate'); % Construct the activate script path
command = sprintf('source %s && python3.11 %s %s', activate_script, python_script, arguments);

[status, result] = system(command);
disp(result)

