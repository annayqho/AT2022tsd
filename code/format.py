# Formatting for plots

def get_format():
    dict = {}

    # Font sizes
    dict['font_large'] = 13
    dict['font_med'] = 11
    dict['font_small'] = 9

    # Colors
    dict['colors'] = {}
    dict['colors']['3'] = ['#003f5c', '#bc5090', '#ffa600']
    dict['colors']['4'] = ['#003f5c', '#7a5195', '#ef5675', '#ffa600']
    dict['colors']['5'] = ['#003f5c', '#58508d', '#bc5090', '#ff6361', '#ffa600']
    dict['colors']['6'] = ['grey','#003f5c', '#58508d', '#bc5090', '#ff6361', '#ffa600']
    dict['colors']['9'] = ['k', '#003f5c', '#2f4b7c', '#665191', '#a05195', '#d45087', '#f95d6a', '#ff7c43', '#ffa600']
    dict['colors']['11'] = ['#ffa600', 'grey', 'k', '#003f5c', '#2f4b7c', '#665191', '#a05195', '#d45087', '#f95d6a', '#ff7c43', '#ffa600']

    return dict
