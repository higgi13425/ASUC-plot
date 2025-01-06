# README -----------------------------------------------------------------------
# 
#   file: analysis-functions.R
#   goal: centralize all functions used for analysis of ASUC cohort
#   
# ------------------------------------------------------------------------------

# packages ---------------------------------------------------------------------

  library(glue)

# univariable analysis ---------------------------------------------------------

  ## calculate ORs for list of continuous predictors 

    get_LR_summary = function(df, var_list, output_var){
      
      output_df = data.frame()
      
      for(var in var_list){
        
        # generate formula object
        temp_formula = as.formula(paste(output_var, ' ~ ', var))
        
        # fit logistic regression model
        temp_mod     = glm(temp_formula, family = binomial(link = 'logit'), data = df)
        
        # gather coefficient info and calculate OR + 95% CI
        temp_summary = data.frame(summary(temp_mod)$coefficients)[2, ] %>%
          mutate(
            lower_ci = Estimate - 1.96 * Std..Error, 
            upper_ci = Estimate + 1.96 * Std..Error, 
            OR = exp(Estimate), 
            lower_OR = exp(lower_ci), 
            upper_OR = exp(upper_ci)
          )
        
        # organize into dataframe and bind to existing data
        temp_coef     = round(temp_summary[1, 1], 4)
        temp_p        = round(temp_summary[1, 4], 7)
        temp_OR       = round(temp_summary[1, 7], 2)
        temp_lower_OR = round(temp_summary[1, 8], 2)
        temp_upper_OR = round(temp_summary[1, 9], 2)
        
        temp_output = data.frame(
          variable    = var, 
          coefficient = temp_coef, 
          p           = temp_p, 
          OR          = temp_OR, 
          lower_OR    = temp_lower_OR, 
          upper_OR    = temp_upper_OR
        )
        
        temp_output = temp_output %>% 
          mutate(
            OR_summary = glue("{OR} ({lower_OR}, {upper_OR})")
          ) %>% 
          select(variable, p, OR_summary)
        
        output_df = rbind(output_df, temp_output)
        
      }
      
      output_df = output_df %>% 
        arrange(p) %>% 
        mutate(
          p = ifelse(
            p < 0.001, "< 0.001", 
            round(p, digits = 3) %>% as.character() 
          )
        )
      
      return(output_df)
      
    }
    
