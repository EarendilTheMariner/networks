function W = dalien_network(size,f_exc,sparsity)
    % generate dalien network of size (size), fraction of excitatory edges
    % (f_exc) and sparsity (sparsity)

    n = size;
    W = zeros(n);
    p = sparsity; 
    
    neuron_types = [ones(1, n*f_exc), -1*ones(1,(n*(1-f_exc)))]; % 1xn array signifying neuron types (1 or -1)
    
    for i = 1:n % loop through rows (outgoing edges)
        for j = 1:n % loop through entries of row i
            if i ~= j && rand <= p % make edge from neuron i to j with probability p
                W(i, j) = neuron_types(i) * rand; % weight with correct sign and between 0 and (+/-)1
            end
        end
    end
  
end

