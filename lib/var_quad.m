function [ VAR_match ] = var_quad( params,setup,store_responses )
%objective function for VAR responses matching
[ MA_matrices ] = wrap_BM_var_opt( params, setup );
MA_h=MA_matrices(:,:,1:11);
% size(MA_h)
% size(store_responses)
VAR_match=(MA_h(:)-store_responses(:))'*(MA_h(:)-store_responses(:));
end

