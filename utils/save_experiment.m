function save_experiment(info,experiment_path,param)
%SAVE_EXPERIMENT Save experiment
%   Save experiment in specified path along with animation of its evolution

% Check if directory exists, create otherwise
if experiment_path(end) ~= '/', experiment_path = [experiment_path '/']; end
if (exist(experiment_path, 'dir') == 0), mkdir(experiment_path); end

if ~isfield(param,'operator')
    param.operator = IdentityOperator(size(info.u_history(:,:,end)));
end

% Save experiment detailed information
save(experiment_path,'info');

% Create a gif animation of image evolution
image_evolution_path = [experiment_path 'image_evolution.gif'];
create_animation_image(info.u_history,image_evolution_path);

% Create a gif animation of the parameter evolution
parameter_evolution_path = [experiment_path 'parameter_evolution.gif'];
create_animation_parameter(info.sol_history,parameter_evolution_path,param.operator);

end

