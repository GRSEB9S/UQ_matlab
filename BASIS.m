%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PsiBasis] =   BASIS(nvars, OrderPol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ->Par's routine (unchanged).

%     Generate array PsiBasis{n} , n=1, PsiBasisSize
%     with PsiBasis{n} = array of NbRva elements summing up to OrderPol
%
% <=> Generate lists t / t[i] = sum_{j=1}^{i} (PsiBasis{n}(i)+1)
%
%   Example NbRva = 4 , OrderPol = 2
%   PsiBasis{n}  : [0 0 0 2]   ->  t : [1 2 3]
%               [0 1 0 1]   ->      [1 3 4]
%               [2 0 0 0]   ->      [3 4 5]
%
%   In the sequel tables t are generated, and converted to PsiBasis{n}



%fprintf('\n  * Computing the polynomial chaos basis Psi_j :');

%%------------------------------------------------------------%%
%      Compute the number of basis vectors Psi
%      in the Polynomial chaos associated with
%      (NbRva , OrderPol)
%%------------------------------------------------------------%%

NbRva = nvars;

PsiBasisSize = 0 ;
for k = 0 : OrderPol
    PsiBasisSize = PsiBasisSize + factorial(NbRva + k -1)/ ...
        (factorial(k) * factorial (NbRva-1));
end
%fprintf(' (size = %3.0d)  ', PsiBasisSize );

%%------------------------------------------------------------%%
%      Generate the basis vectors
%%------------------------------------------------------------%%
%PsiBasis = cell(PsiBasisSize,1);  % Basis numbering starting with 1
PsiBasis = zeros(int16(PsiBasisSize),NbRva);
%PsiBasis{1} = zeros([1 NbRva]);

% Initialize the generation
MM = NbRva - 1 ;
n = 1;

% Particular case
if (NbRva > 0)
    switch NbRva
        
        case 1 % One-dimensional Polynomial Chaos <=> Hermite Polynomials
            for i = 1 : OrderPol
                %PsiBasis{i+1}(1) = i;
                PsiBasis(i+1,1) = i;
            end;
            
        otherwise % Higher dimensional Polynomial Chaoses
            for CurrentOrder = 1 : OrderPol
                EndGenere = 0;
                FirstThisOrder = 0;
                
                while (EndGenere == 0)
                    n = n +1 ;
                    % First list t for order CurrentOrder
                    if ( FirstThisOrder == 0)
                        for i=1 : MM
                            t(i) = i;
                        end;
                        FirstThisOrder =1 ;
                    else
                        % Regular incrementation
                        if (t(MM) < (MM + CurrentOrder))
                            t(MM) = t(MM) + 1;
                        else  % t(MM) = tmax = MM +CurrentOrder
                            j = MM;
                            while (t(j) == j + CurrentOrder )
                                j = j - 1 ;
                            end;
                            t(j) = t(j) + 1 ;
                            for k =(j + 1) :  MM
                                t(k) = t(j) + k - j ;
                            end;
                        end;
                    end;
                    
                    % Direct Translating t into PsiBasis{n}
                    %PsiBasis{n}(1) = t(1) -1;
                    PsiBasis(n,1) = t(1) -1;
                    for i=2 : MM
                        %PsiBasis{n}(i) = t(i) - t(i-1) -1 ;
                        PsiBasis(n,i) = t(i) - t(i-1) -1;
                    end;
                    %PsiBasis{n}(NbRva) = NbRva + CurrentOrder - t(MM) -1 ;
                    PsiBasis(n,NbRva) = NbRva + CurrentOrder - t(MM) -1;
                    
                    % End of generation of order CurrentOrder
                    if (t(1) == (CurrentOrder+1))
                        EndGenere = EndGenere + 1;
                    end;
                end;
            end;
    end;
else
    disp('The number of random variable has to be greater than 0');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
