

function [sym_fname] = write_dx_sym(E, Pm, M, Y, D, RefNode, ode_fname)
[Sl_star, Vl, genbus_id, loadbus_id, xd, V_0, V_0_A, Sg_star, NOB, D, H, M, Gn, ns, T_end, Ts]=indi_conds
AllNodes = [1:Gn];
NonRefNodes = setdiff(AllNodes,RefNode);

G = real(Y);
B = imag(Y);

%==================================================
% Create a blank m-file to write the ODE function
%==================================================

sym_fname = strcat(ode_fname,'_sym');

fid = fopen(strcat(sym_fname,'.m'),'w'); % Creates a new file with the SYMBOLIC dynamics

%=============================================================
% WRITE in the newly created m-file with the file handle 'fid'
%=============================================================

% First few lines
fprintf(fid, 'function [f,x] = dx_for_9bus_sym \n')

fprintf(fid, '\n%%This funtion creates the function handle describing the SYM dynamics of 9-bus system.\n')
define_symvar = 'syms ';

for p = 1:length(NonRefNodes)
    define_symvar = [define_symvar sprintf('x%d x%d ',p,p+length(NonRefNodes))];
end
fprintf(fid, '\n\n%s real\n',define_symvar)


% Then the dx(1) to dx(2)

for p = 1:length(NonRefNodes)
    fprintf(fid, '\nf%d = x%d;', p, p+length(NonRefNodes));
end

%===============================
% Finally, the dx(3) to dx(4)
%===============================

% First, compute the 'strings' sum_str{1} to sum_str{2}, and dd3_str
%The following summation will calculate the expression for Pe=electric
%plower, as in the swing ew as ddelta/dt = Pm - Pe - Pd. <=====This Pe for
%each generator i=1 to 2. The same for 3 is calculated seperatety...
%Reason::: We have not included delta3 and omega3 in our original state
%vector. So it is little difficult to put in the same loop.
%------the swing equtions - similar to this one
for p = 1:length(NonRefNodes)
    i = NonRefNodes(p);
    sum_str{i} = sprintf(' ( %.6g * ( %.6g * sin(x%d-0) + %.6g * cos(x%d-0))', E(i,1)*E(RefNode,1), B(i,RefNode), p, G(i,RefNode), p); %We use the term corresponding to delta_3 as initial value for sum_i.
    for q = 1:length(NonRefNodes)
        j = NonRefNodes(q);
        sum_str{i} = strcat(sum_str{i}, sprintf(' + %.6g * (%.6g * sin(x%d-x%d) + %.6g * cos(x%d-x%d))', E(i,1)*E(j,1), B(i,j), p, q, G(i,j), p, q)); % instead of this one we need actual equation.
    end
    sum_str{i} = strcat(sum_str{i},')');
end

%Now let us find the dynamics corresponding to delta3 an omega3.
sum_str{RefNode} = sprintf(' ( %.6g', E(RefNode,1)^2*G(RefNode,RefNode)); %sum for ref RefNode
for p = 1:length(NonRefNodes)
    i = NonRefNodes(p);
    sum_str{RefNode} = strcat(sum_str{RefNode}, sprintf(' + %.6g * (%.6g * sin(0 - x%d) + %.6g * cos(0 - x%d))', E(RefNode,1)*E(i,1), B(RefNode,i), p, G(RefNode,i), p));
end
sum_str{RefNode} = strcat(sum_str{RefNode},')');

dd3_str = sprintf('(1/ %.6g) .* (%.6g - %s)', M(RefNode), Pm(RefNode,1), sum_str{RefNode}); %dd3 is the reference i.e., (dw3/dt) when, w3=0,  delta3=0.

% Second, write these strings into the dx-equations
for p = 1:length(NonRefNodes)
    i = NonRefNodes(p);
    dx_str = sprintf('(1/ %.6g) .* (%.6g - %s - %.6g * x%d) - %s', M(i), Pm(i,1), sum_str{i}, D(i,1), (p+length(NonRefNodes)), dd3_str);
    fprintf(fid, '\nf%d = %s;', (p+length(NonRefNodes)), dx_str);
end

define_fsym = 'f = [';
define_xsym = 'x = [';
for p = 1:2*length(NonRefNodes)
    define_fsym = [define_fsym sprintf('f%d; ',p)];
    define_xsym = [define_xsym sprintf('x%d; ',p)];
end
define_fsym = strcat(define_fsym, '];');
define_xsym = strcat(define_xsym, '];');
fprintf(fid, '\n\n%s',define_fsym);
fprintf(fid, '\n\n%s',define_xsym);

fclose(fid);
