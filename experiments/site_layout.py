"""Rearrange the existing procedure blocks without changing their descriptions."""
from bs4 import BeautifulSoup

DETECTION=['Global sum','Maximum entry','Column pooling','Exhaustive scan','Combined sum and scan']
RECOVERY=['Exhaustive block maximization','Spectral projection','Thresholded spectral projection']


def procedure_content(section):
    return {a.h3.get_text(' ',strip=True):dict(text=a.get_text(' ',strip=True),
                formulas=[m['data-tex'] for m in a.select('[data-tex]')])
            for a in section.select('article.method-card')}


def reorganize_methods(source):
    soup=BeautifulSoup(source,'html.parser')
    section=soup.find('section',id='methods')
    before=procedure_content(section)
    container=section.find('div',class_='container')
    heading=container.find('div',class_='section-head').extract()
    branch=container.find('div',class_='inference-branch').extract()
    leaves=branch.select('.branch-leaves strong')
    leaves[0]['id']='detection-branch-title';leaves[1]['id']='recovery-branch-title'
    cards={a.h3.get_text(' ',strip=True):a.extract() for a in section.select('article.method-card')}
    recovery_intro=container.find('div',id='recovery-procedures').find('p',recursive=False).extract()
    boundary=container.find('div',class_='recovery-boundary').extract()
    comparison=container.find('div',id='procedure-comparison').extract()
    # This directional reference changes only because the equation now precedes the table.
    for text in comparison.find_all(string=True):
        if 'displayed below.' in text:text.replace_with(text.replace('displayed below.','given above.'))
    leaves=branch.find('div',class_='branch-leaves').extract()
    container.clear();container.append(heading)
    flow=soup.new_tag('div',attrs={'class':'inference-flow'})
    flow.append(branch);flow.append(leaves)
    tree=soup.new_tag('div',attrs={'class':'procedure-tree','aria-label':'Detection and recovery procedure branches'})
    for goal,names in [('detection',DETECTION),('recovery',RECOVERY)]:
        lane=soup.new_tag('div',attrs={'class':'procedure-lane '+goal+'-lane',
              'aria-labelledby':goal+'-branch-title','id':goal+'-procedures'})
        mobile_heading=soup.new_tag('h3',attrs={'class':'lane-mobile-title'})
        mobile_heading.string=goal.capitalize();lane.append(mobile_heading)
        if goal=='recovery':
            recovery_intro['class']=['lane-intro'];lane.append(recovery_intro)
        for number,name in enumerate(names,1):
            node=soup.new_tag('div',attrs={'class':'procedure-node','data-step':f'{number:02d}'})
            node.append(cards[name]);lane.append(node)
        end=soup.new_tag('div',attrs={'class':'lane-end'})
        end.string='Establish presence' if goal=='detection' else 'Recover both index sets'
        lane.append(end);tree.append(lane)
    flow.append(tree);container.append(flow)
    title=soup.new_tag('h3',attrs={'id':'recovery-boundary-title','class':'recovery-boundary-title'})
    title.string='The signal needed for exact recovery'
    boundary.insert(0,title);boundary['aria-labelledby']='recovery-boundary-title'
    container.append(boundary)
    title=soup.new_tag('h3',attrs={'class':'comparison-title'})
    title.string='Compare the procedures';comparison.insert(0,title);container.append(comparison)
    assert before==procedure_content(section),'A procedure description or formula changed.'
    return str(soup)
