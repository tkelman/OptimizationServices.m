clear
instance = OSinstance;
instance.instanceHeader.name = 'Bonmin Example';
instance.instanceHeader.source = 'Bonmin example folder';
instance.instanceHeader.description = sprintf(['Objective value: -1.70711\n' ...
    '\t\t\tSolution:\n' ...
    '\t\t\tx[0] = 0, x[1] = 0.853553, x[2] = 0.853553, x[3] = 1']);

% var = instance.createVariable(name, type, lb, ub, mult)
x0 = instance.createVariable('x0', 'B', 0, 1);
x1 = instance.createVariable('x1', 'C', 0);
x2 = instance.createVariable('x2', 'C', 0);
x3 = instance.createVariable('x3', 'I', 0, 5);

% instance.setObjective(expression, name, maxOrMin)
instance.setObjective(x0 - x1 - x2, 'minCost', 'min');

% instance.createConstraint(expression, name, constant, lb, ub, mult)
instance.createConstraint((x1 - 0.5)^2 + (x2 - 0.5)^2, [], [], [], 0.25);
instance.createConstraint(x0 - x1,                     [], [], [], 0);
instance.createConstraint(x1 + x2 + x3,                [], [], [], 2);

xmlwrite(instance.document)
xmlwrite('bonminEx1_nonlinear_out.osil', instance.document);
