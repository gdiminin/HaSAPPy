from scipy.stats  import fisher_exact

def main(GroupAnalysis,DATA):
    for group in GroupAnalysis.Others.name:
        name_control = '%s_KI_sum' %GroupAnalysis.Reference.name
        name_selected = '%s_KI_sum' %group
        total_KI_control = DATA[name_control].sum() - DATA.ix['_no_feature',name_control]
        total_KI_selected = DATA[name_selected].sum() - DATA.ix['_no_feature',name_selected]
        temporary_dataframe = DATA.copy()
        temporary_dataframe.ix['_no_feature',name_control] =1
        temporary_dataframe.ix['_no_feature',name_selected] =1

        DATA[‘%s_Score_fisher’% group] =  temporary_dataframe.apply(lambda r: fisher_exact([[r[name_selected], r[name_control]],[(total_KI_selected-r[name_selected]),(total_KI_control-r[name_control])]],alternative = 'greater')[1],axis=1)
        DATA.ix['_no_feature’,’%s_Score_fisher’% group] = 1
    return DATA
