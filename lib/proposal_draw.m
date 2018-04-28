function [ draw ] = proposal_draw( posterior_draw,cov_matrix, scaling,setup )
%gives draw from a RW Gaussian propsal density with mean posterior_draw and
%covariance matrix scaling*cov_matrix

draw=mvnrnd(posterior_draw,scaling*cov_matrix);



end

