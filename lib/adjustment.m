function [ adjust_factor ] = adjustment( old_draw,proposal,setup )
%This function adjusts the MH acceptance probability for the fact that some proposal processes for parameters are transformed so that the parameters are restricted


adjust_factor=0; %in logs

if setup.proposal==1 && ((setup.length_log+setup.length_logit+setup.length_logit_general)>0) %for now this only works for the standard RW MH algorithm
    
    if setup.length_log>0
        for jj=1:setup.length_log
            adjust_factor=adjust_factor+proposal(setup.index_log(jj))-old_draw(setup.index_log(jj));
            
        end
        
        
    end
    
    if (setup.length_logit+setup.length_logit_general)>0
    old_draw_rescaled = inv_transform( old_draw,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );
    proposal_rescaled = inv_transform( proposal,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );

        
    end
    
    
    
    
     if setup.length_logit>0
        for jj=1:setup.length_logit
            adjust_factor=adjust_factor+log(proposal_rescaled(setup.index_logit(jj))-proposal_rescaled(setup.index_logit(jj))^2)-log(old_draw_rescaled(setup.index_logit(jj))-old_draw_rescaled(setup.index_logit(jj))^2);
            
        end
        
        
     end
    
    
      if setup.length_logit_general>0
        for jj=1:setup.length_logit_general
            p1=(proposal_rescaled(setup.index_logit_general(jj))-setup.logit_general_lb(jj))/(setup.logit_general_ub(jj)-setup.logit_general_lb(jj));
            p2=(old_draw_rescaled(setup.index_logit_general(jj))-setup.logit_general_lb(jj))/(setup.logit_general_ub(jj)-setup.logit_general_lb(jj));
  
            adjust_factor=adjust_factor+log(p1-p1^2)-log(p2-p2^2);
            
        end
        
        
    end
     
     
    
end


end

