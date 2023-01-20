# Concatenate, sort, etc
radio = pd.concat([vla,noema],axis=0,ignore_index=True).sort_values('Date', ignore_index=True)
