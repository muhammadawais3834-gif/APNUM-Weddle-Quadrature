%% APNUM_Appendix_NumericalCode.m
% Reproducible numerical experiments for the APNUM manuscript.
% Methods:
%   NC-2 : composite Simpson (closed Newton-Cotes) rule
%   G-2  : composite two-point Gauss-Legendre rule
%   M-2  : (2/5) NC-2 + (3/5) G-2
%
% Test problems:
%   E1: f(x)=x^6,                         [0,1]
%   E2: f(x)=e^(4x),                      [0,1]
%   E3: f(x)=-sin(x),                     [0,pi/2]
%   E4: f(x)=-(1/pi^2)sin(pi*x),          [0,1]
%   E5: q-digamma psi_q(x), q=0.5,        [1,2]
%   E6: f(x)=W'_nu(x), normalized Bessel, [1,2]
%
% MATLAB R2016b or later.

clear; clc; close all;
format long e;

%% Reproducibility parameters
Nlist = [12 24 48 96 192 384];
qpar  = 0.5;

% IMPORTANT: this value must match the final manuscript everywhere.
% The original MATLAB source supplied for the manuscript uses nu = 0.5.
nu = 0.5;

plotFloor = 1e-16;

%% Define examples
examples = struct([]);

examples(1).id = 'E1';
examples(1).name = 'f(x)=x^6';
examples(1).f = @(x) x.^6;
examples(1).a = 0;
examples(1).b = 1;
examples(1).Iref = 1/7;
examples(1).titleTex = '$f(x)=x^6$';

examples(2).id = 'E2';
examples(2).name = 'f(x)=e^(4x)';
examples(2).f = @(x) exp(4*x);
examples(2).a = 0;
examples(2).b = 1;
examples(2).Iref = (exp(4)-1)/4;
examples(2).titleTex = '$f(x)=e^{4x}$';

examples(3).id = 'E3';
examples(3).name = 'f(x)=-sin(x)';
examples(3).f = @(x) -sin(x);
examples(3).a = 0;
examples(3).b = pi/2;
examples(3).Iref = -1;
examples(3).titleTex = '$f(x)=-\sin(x)$';

examples(4).id = 'E4';
examples(4).name = 'f(x)=-(1/pi^2)sin(pi*x)';
examples(4).f = @(x) -(1/pi^2)*sin(pi*x);
examples(4).a = 0;
examples(4).b = 1;
examples(4).Iref = -2/pi^3;
examples(4).titleTex = '$f(x)=-\frac{1}{\pi^2}\sin(\pi x)$';

examples(5).id = 'E5';
examples(5).name = sprintf('q-digamma, q=%.3g',qpar);
examples(5).f = @(x) qdigamma_apnum(x,qpar);
examples(5).a = 1;
examples(5).b = 2;
examples(5).Iref = 0;
examples(5).titleTex = sprintf('$q$-digamma, $q=%.1f$',qpar);

examples(6).id = 'E6';
examples(6).name = sprintf('f(x)=W''_nu(x), nu=%.3g',nu);
examples(6).f = @(x) W1_apnum(x,nu);
examples(6).a = 1;
examples(6).b = 2;
examples(6).Iref = W0_apnum(2,nu)-W0_apnum(1,nu);
examples(6).titleTex = sprintf('$f(x)=\\mathcal{W}_{\\nu}''(x)$, $\\nu=%.1f$',nu);

%% Compute approximations and errors
methodNames = {'NC-2','G-2','M-2'};
nExamples = numel(examples);
nN = numel(Nlist);
nMethods = numel(methodNames);

Approx = zeros(nExamples,nN,nMethods);
AbsErr = zeros(nExamples,nN,nMethods);
RelErr = NaN(nExamples,nN,nMethods);
ObsOrd = NaN(nExamples,nN,nMethods);

for e = 1:nExamples
    f = examples(e).f;
    a = examples(e).a;
    b = examples(e).b;
    Iref = examples(e).Iref;

    fprintf('\n%s : %s\n',examples(e).id,examples(e).name);
    fprintf('Reference integral = %.15e\n',Iref);

    for j = 1:nN
        N = Nlist(j);

        Qnc = composite_NC2(f,a,b,N);
        Qg  = composite_G2(f,a,b,N);
        Qm  = (2/5)*Qnc + (3/5)*Qg;

        Approx(e,j,1) = Qnc;
        Approx(e,j,2) = Qg;
        Approx(e,j,3) = Qm;

        for m = 1:nMethods
            AbsErr(e,j,m) = abs(Approx(e,j,m)-Iref);
            if Iref ~= 0
                RelErr(e,j,m) = AbsErr(e,j,m)/abs(Iref);
            end
        end
    end

    for m = 1:nMethods
        for j = 1:nN-1
            if AbsErr(e,j,m) > 0 && AbsErr(e,j+1,m) > 0
                ObsOrd(e,j,m) = log2(AbsErr(e,j,m)/AbsErr(e,j+1,m));
            end
        end
    end

    fprintf('      N          NC-2 error          G-2 error           M-2 error\n');
    for j = 1:nN
        fprintf('%7d    % .12e    % .12e    % .12e\n', ...
            Nlist(j),AbsErr(e,j,1),AbsErr(e,j,2),AbsErr(e,j,3));
    end
end

%% Export complete numerical data
rows = nExamples*nN*nMethods;
Example = strings(rows,1);
Function = strings(rows,1);
Ncolumn = zeros(rows,1);
Method = strings(rows,1);
Reference = zeros(rows,1);
Approximation = zeros(rows,1);
AbsoluteError = zeros(rows,1);
RelativeError = NaN(rows,1);
ObservedOrder = NaN(rows,1);

r = 0;
for e = 1:nExamples
    for j = 1:nN
        for m = 1:nMethods
            r = r+1;
            Example(r) = string(examples(e).id);
            Function(r) = string(examples(e).name);
            Ncolumn(r) = Nlist(j);
            Method(r) = string(methodNames{m});
            Reference(r) = examples(e).Iref;
            Approximation(r) = Approx(e,j,m);
            AbsoluteError(r) = AbsErr(e,j,m);
            RelativeError(r) = RelErr(e,j,m);
            ObservedOrder(r) = ObsOrd(e,j,m);
        end
    end
end

Results = table(Example,Function,Ncolumn,Method,Reference,Approximation, ...
    AbsoluteError,RelativeError,ObservedOrder);
writetable(Results,'APNUM_quadrature_results.csv');

%% Generate the six 2D convergence plots
for e = 1:nExamples
    figure('Color','w','Position',[100 100 900 580]);

    yNC = max(squeeze(AbsErr(e,:,1)),plotFloor);
    yG  = max(squeeze(AbsErr(e,:,2)),plotFloor);
    yM  = max(squeeze(AbsErr(e,:,3)),plotFloor);

    loglog(Nlist,yNC,'-o','LineWidth',1.6,'MarkerSize',6); hold on;
    loglog(Nlist,yG,'--s','LineWidth',1.6,'MarkerSize',6);
    loglog(Nlist,yM,'-.^','LineWidth',1.6,'MarkerSize',6);

    grid on; box on;
    xlabel('Nominal resolution $N$','Interpreter','latex');
    ylabel('Absolute quadrature error','Interpreter','latex');
    title(examples(e).titleTex,'Interpreter','latex');

    legend({'Newton--Cotes NC-2 (Simpson)', ...
            'Gauss--Legendre G-2', ...
            'Mixed M-2'}, ...
            'Location','best','Interpreter','latex');

    xticks(Nlist);
    xticklabels(string(Nlist));

    pdfName = sprintf('%s_verified_2D_comparison.pdf',examples(e).id);
    pngName = sprintf('%s_verified_2D_comparison.png',examples(e).id);

    exportgraphics(gcf,pdfName,'ContentType','vector');
    exportgraphics(gcf,pngName,'Resolution',600);
end

%% Print compact manuscript values
fprintf('\nCompact manuscript values: N=12, N=24 and p_obs\n');
for e = 1:nExamples
    fprintf('\n%s\n',examples(e).id);
    for m = 1:nMethods
        fprintf('%s: E12 = %.12e, E24 = %.12e, p = %.6f\n', ...
            methodNames{m},AbsErr(e,1,m),AbsErr(e,2,m),ObsOrd(e,1,m));
    end
end

%% Local functions
function Q = composite_NC2(f,a,b,N)
% Composite Simpson (NC-2) rule.
    if mod(N,2) ~= 0
        error('NC-2 requires an even N.');
    end

    nPanels = N/2;
    H = (b-a)/nPanels;
    Q = 0;

    for k = 0:nPanels-1
        left = a+k*H;
        right = left+H;
        mid = (left+right)/2;

        Q = Q + (right-left)/6 * ...
            (f(left)+4*f(mid)+f(right));
    end
end

function Q = composite_G2(f,a,b,N)
% Composite two-point Gauss-Legendre rule.
    if mod(N,2) ~= 0
        error('G-2 comparison requires an even N.');
    end

    nPanels = N/2;
    H = (b-a)/nPanels;
    xi = 1/sqrt(3);
    Q = 0;

    for k = 0:nPanels-1
        left = a+k*H;
        right = left+H;
        mid = (left+right)/2;
        half = (right-left)/2;

        x1 = mid-half*xi;
        x2 = mid+half*xi;

        Q = Q + half*(f(x1)+f(x2));
    end
end

function y = qdigamma_apnum(x,q)
% q-digamma for 0<q<1:
% psi_q(x)=-log(1-q)+log(q)*sum q^(k+x)/(1-q^(k+x)).
    if ~(q > 0 && q < 1)
        error('This implementation assumes 0 < q < 1.');
    end

    y = zeros(size(x));
    tol = 1e-15;
    maxK = 100000;

    for j = 1:numel(x)
        s = 0;
        for k = 0:maxK
            r = q^(k+x(j));
            term = r/(1-r);
            s = s+term;
            if abs(term) < tol
                break;
            end
        end

        if k == maxK
            warning('q-digamma series reached maxK.');
        end

        y(j) = -log(1-q)+log(q)*s;
    end
end

function W = W0_apnum(x,nu)
% Normalized modified Bessel function:
% W_nu(x)=2^nu*Gamma(nu+1)*x^(-nu)*I_nu(x).
    W = 2^nu*gamma(nu+1).*x.^(-nu).*besseli(nu,x);
end

function Wp = W1_apnum(x,nu)
% W'_nu(x)=x*W_(nu+1)(x)/(2*(nu+1)).
    Wp = x.*W0_apnum(x,nu+1)./(2*(nu+1));
end
