function jordan_matrix = jordan(A, tolerance)

if nargin < 2, tolerance = 0.005; end

[eigenvalues,mult] = heltalsev(A,tolerance)

% Eliminates duplicates
ev_unique = unique(round(eigenvalues));

n_rows = size(A,1);

jordan_matrix = zeros(n_rows);
n_unique_ev = size(ev_unique,1);
eye_matrix = eye(n_rows);
position = 1;               %Diagonal position of a11 in inserting block.

for ii = 1:size(mult,1)

    if mult(ii) == 1
        jordan_matrix(position, position) = ev_unique(ii);
        position = position + 1;
        
    elseif mult(ii) == 2
        
        B = A - eye_matrix*ev_unique(ii);
        if abs(B) < tolerance
            B = round(B);
        end
        
        r1 = rank(B);
        p1 = n_rows - r1;
        
        if p1 == 2
            jordan_matrix(position, position) = ev_unique(ii);
            jordan_matrix(position + 1, position + 1) = ev_unique(ii);
            position = position + 2;
        else
            D = [ev_unique(ii) 1; 0 ev_unique(ii)];
            jordan_matrix(position: position + 1, position: position + 1) = D;
            position = position + 2;
        end

    else
        % every k calculates one jordan table for unique eigenvalue
        for k = 1:n_unique_ev 
            
            % terminate when k reaches multiplicity or b = 0
            k_not_multiplicity = true;          
            
            index = 1;
            while(k_not_multiplicity)
                
                B = A - eye_matrix*eigenvalues(k);
                B = B^index;
          
                % Tol works for 0,1,3-6 testmatris
                tolerance = norm(B)*10;   
                
                if abs(B) < tolerance
                    B = round(B);
                end
                %Rank of (A-lambda*I)^i
                rank_matrix(index,k) = rank(B);     
                
                %p_matrix is dim(Ker(A-lambda*I(^i))
                p_matrix(index,k) = n_rows - rank_matrix(index,k);              
                if index == 1
                    %initiation of b1
                    b(index,k) = p_matrix(index,k);                        
                else
                    % b is defined like this, b defines stop condition for the loop
                    b(index,k) = p_matrix(index,k) - p_matrix(index-1,k);
                end
                if b(index,k) == 0
                    % Stop condition
                    k_not_multiplicity = false;                                
                end
                if index == mult(k)
                    % Stop condition
                    k_not_multiplicity = false;                               
                end
                index = index + 1;
            end
        end
        form_jordan_block = zeros(size(b));
        b(1:end,:);
        
        % n determines form of Jordan blocks
        form_jordan_block(1:end-1,:) = b(1:end-1,:) - b(2:end,:);  
        form_jordan_block(end,:) = b(end,:);
        
        % for loop inserting blocks
        for k = 1:size(form_jordan_block,2)                             
            for index = 1:size(form_jordan_block,1)
                if form_jordan_block(index,k) > 0
                    for antal = 1:form_jordan_block(index,k)
                        C = ones(index);
                        
                        %superdiagonal of 1´s if needed
                        D = triu(C,1) - triu(C,2);
                  
                        D = D + eye(index) * ev_unique(k);
                        jordan_matrix(position:position+size(D,1)-1,position:position+size(D,1)-1) = D;
                        position = position + index;
                    end
                end
            end
        end
    end
end
end
