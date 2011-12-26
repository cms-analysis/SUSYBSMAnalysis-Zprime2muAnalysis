import copy
from FWCore.PythonUtilities.LumiList import LumiList

def for_cmssw(ll):
    return ll.getCMSSWString().split(',')

# These numbers dictate how the rereco, prompt, DCS-only jsons are
# combined below.
first_run = 160404
last_rereco_run = 163869
last_prompt_run = 180252
last_run = 180252
assert last_prompt_run > last_rereco_run
assert last_run > first_run

# Sometimes the same run-range json gets made in other versions.
prompt_version = ''

# Lumis to manually throw out.
# This first list amounts to 12/pb of data in which the prescale changed mid-LS. See https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/1558.html
to_remove = { "162811": [[4, 4]], "162909": [[3, 3]], "162929": [[131, 131], [273, 273]], "163233": [[274, 274]], "163252": [[4, 4]], "163286": [[108, 108]], "163296": [[3, 3]], "163332": [[7, 7]], "163374": [[671, 671]], "163589": [[46, 46]], "163759": [[182, 182]], "163795": [[24, 24], [33, 33]], "163817": [[118, 118]], "165098": [[78, 78]], "165120": [[27, 27]], "165121": [[242, 242]], "165364": [[1132, 1132]], "165415": [[812, 812]], "165467": [[113, 113]], "165472": [[493, 493]], "165506": [[57, 57]], "165514": [[40, 40]], "165523": [[5, 5]], "165537": [[13, 13]], "165542": [[161, 161]], "165567": [[540, 540]], "165570": [[696, 696]], "165633": [[2, 2]], "165970": [[9, 9]], "165993": [[773, 773]], "166033": [[13, 13], [27, 27]], "166034": [[203, 203]], "166049": [[15, 15], [57, 57]], "166346": [[7, 7]], "166374": [[3, 3]], "166380": [[637, 637]], "166408": [[29, 29]], "166429": [[4, 4]], "166438": [[1, 1]], "166462": [[15, 15]], "166502": [[2, 2]], "166514": [[12, 12]], "166530": [[37, 37], [42, 42]], "166565": [[12, 12]], "166701": [[79, 79]], "166781": [[2, 2]], "166787": [[12, 12]], "166839": [[5, 5]], "166842": [[7, 7]], "166859": [[26, 26]], "166864": [[318, 318]], "166895": [[135, 135]], "166960": [[60, 60]], "167098": [[423, 423], [447, 447]], "167102": [[6, 6]], "167282": [[38, 38]], "167284": [[891, 891]], "167673": [[2, 2]], "167675": [[517, 517]], "167740": [[4, 4]], "167807": [[1406, 1406]], "167830": [[171, 171], [458, 458]], "167898": [[4, 4], [1334, 1334], [1655, 1655]], "167913": [[12, 12], [14, 14], [360, 360]], "169811": [[11, 11]], "170040": [[71, 71], [112, 112]], "170249": [[41, 41], [45, 45]], "170286": [[3, 3], [83, 83], [87, 87]], "170348": [[6, 6], [105, 105], [113, 113], [143, 143]], "170354": [[290, 290]], "170376": [[53, 53], [56, 56]], "170382": [[8, 8], [10, 10]], "170452": [[22, 22], [75, 75], [79, 79]], "170527": [[13, 13]], "170722": [[116, 116], [119, 119], [219, 219]], "170759": [[14, 14], [336, 336]], "170876": [[180, 180], [415, 415]], "171050": [[53, 53], [92, 92]], "171106": [[29, 29]], "171156": [[211, 211]], "171178": [[149, 149], [731, 731]], "171369": [[61, 61]], "171446": [[3, 3], [375, 375]], "171484": [[79, 79], [358, 358]], "171578": [[347, 347], [696, 696]], "171812": [[323, 323]], "171876": [[382, 382]], "171897": [[295, 295]], "171921": [[4, 4]], "172033": [[256, 256]], "172163": [[543, 543]], "172268": [[2, 2]], "172389": [[4, 4]], "172399": [[1, 1]], "172411": [[3, 3]], "172485": [[2, 2]], "172488": [[13, 13]], "172495": [[115, 115]], "172630": [[11, 11]], "172778": [[2, 2]], "172791": [[67, 67], [664, 664]], "172799": [[203, 203]], "172802": [[304, 304]], "172819": [[42, 42]], "172822": [[616, 616], [1922, 1922]], "172847": [[1, 1]], "172865": [[2, 2]], "172868": [[553, 553], [1755, 1755]], "172949": [[2, 2], [933, 933]], "172952": [[680, 680]], "172998": [[49, 49], [65, 65], [80, 80]], "173198": [[665, 665]], "173236": [[4, 4], [75, 75]], "173241": [[16, 16], [759, 759]], "173243": [[13, 13]], "173380": [[4, 4], [199, 199]], "173389": [[283, 283]], "173406": [[22, 22]], "173439": [[307, 307]], "173657": [[4, 4], [56, 56], [68, 68], [90, 90]], "173692": [[2, 2], [89, 89], [927, 927], [2250, 2250]], "175832": [[35, 35], [81, 81]], "175857": [[25, 25], [30, 30]], "175872": [[65, 65], [67, 67]], "175873": [[24, 24]], "175874": [[13, 13], [33, 33], [51, 51]], "175906": [[49, 49], [51, 51], [70, 70]], "175921": [[33, 33], [47, 47], [56, 56], [176, 176], [182, 182]], "175971": [[48, 48], [61, 61]], "175975": [[179, 179]], "175990": [[3, 3], [26, 26], [56, 56]], "176023": [[2, 2], [60, 60], [65, 65]], "176201": [[2, 2], [40, 40], [135, 135]], "176207": [[94, 94]], "176286": [[62, 62], [166, 166]], "176304": [[66, 66], [186, 186]], "176309": [[223, 223], [484, 484], [696, 696], [985, 985]], "176461": [[107, 107], [109, 109]], "176464": [[9, 9]], "176467": [[112, 112]], "176547": [[28, 28]], "176548": [[121, 121], [624, 624]], "176701": [[8, 8]], "176702": [[124, 124], [399, 399]], "176765": [[67, 67], [110, 110]], "176771": [[64, 64], [78, 78], [148, 148]], "176795": [[3, 3], [42, 42]], "176796": [[10, 10], [45, 45]], "176797": [[197, 197]], "176799": [[153, 153]], "176841": [[122, 122]], "176844": [[61, 61]], "176886": [[74, 74], [378, 378]], "176928": [[7, 7]], "176929": [[171, 171]], "176933": [[82, 82]], "176982": [[4, 4]], "177053": [[42, 42], [264, 264], [530, 530]], "177074": [[26, 26], [200, 200], [488, 488]], "177088": [[2, 2]], "177095": [[3, 3]], "177096": [[129, 129]], "177131": [[38, 38]], "177138": [[71, 71]], "177139": [[285, 285], [531, 531]], "177183": [[9, 9], [98, 98]], "177201": [[63, 63], [258, 258], [440, 440]], "177222": [[41, 41], [59, 59]], "177313": [[4, 4], [62, 62]], "177318": [[70, 70], [369, 369]], "177449": [[1, 1], [44, 44], [235, 235]], "177452": [[156, 156]], "177507": [[27, 27], [42, 42], [64, 64], [67, 67]], "177509": [[8, 8], [110, 110]], "177515": [[5, 5]], "177718": [[44, 44], [120, 120], [491, 491]], "177730": [[35, 35], [50, 50], [55, 55], [213, 213], [430, 430]], "177776": [[8, 8], [51, 51]], "177782": [[73, 73], [79, 79], [160, 160]], "177783": [[242, 242]], "177786": [[43, 43]], "177788": [[22, 22], [30, 30]], "177875": [[8, 8], [148, 148], [481, 481]], "178003": [[4, 4], [6, 6]], "178004": [[6, 6], [14, 14]], "178098": [[183, 183], [491, 491]], "178110": [[2, 2], [259, 259]], "178116": [[180, 180]], "178160": [[209, 209]], "178162": [[33, 33]], "178208": [[124, 124]], "178365": [[209, 209], [609, 609]], "178420": [[4, 4], [47, 47], [53, 53]], "178421": [[34, 34], [42, 42], [111, 111], [230, 230], [510, 510]], "178479": [[3, 3], [68, 68], [189, 189], [637, 637]], "178667": [[44, 44]], "178675": [[11, 11]], "178703": [[3, 3], [31, 31], [177, 177]], "178712": [[8, 8]], "178786": [[20, 20], [67, 67], [194, 194]], "178803": [[2, 2], [55, 55], [152, 152]], "178840": [[2, 2], [47, 47], [199, 199]], "178854": [[63, 63]], "178920": [[35, 35], [207, 207]], "178970": [[157, 157]], "178985": [[6, 6]], "179411": [[6, 6], [10, 10], [23, 23], [96, 96]], "179431": [[3, 3]], "179434": [[27, 27], [83, 83], [571, 571]], "179497": [[4, 4], [192, 192]], "179547": [[4, 4], [233, 233]], "179563": [[346, 346]], "179816": [[6, 6], [9, 9], [12, 12], [16, 16], [20, 20], [22, 22], [25, 25]], "179828": [[9, 9], [251, 251]], "179889": [[18, 18], [169, 169]], "179959": [[62, 62]], "180072": [[1, 1], [61, 61], [80, 80]], "180076": [[201, 201], [229, 229]], "180093": [[98, 98], [326, 326]], "180241": [[3, 3], [65, 65], [99, 99], [129, 129]], "180249": [[2, 2]], "180250": [[237, 237]]}
to_remove.update({
    '166530': [[1,105]], # During physics-declared, this run had Mu30_v3 prescaled by 20. https://cmswbm.web.cern.ch/cmswbm/cmsdb/servlet/PrescaleChanges?RUN=166530   https://cmswbm.web.cern.ch/cmswbm/cmsdb/servlet/LumiSections?RUN=166530
    '167102': [[1,7]],   # Mu30_v5 was prescaled by 35 then went down to prescale 1 after these lumis.  https://cmswbm2.web.cern.ch/cmswbm2/cmsdb/servlet/PrescaleChanges?RUN=167102
    })
to_remove = LumiList(compactList=to_remove)

# These runs are <= last_prompt_run, but they were not actually
# considered in the certification for the latest prompt JSON. So,
# don't drop them from the DCS-only list when combining later.
holes = []

runs_to_remove_from_dcsonly = range(first_run, last_prompt_run+1)
for hole in holes:
    print 'goodlumis warning: re-adding "hole" run %i from DCS-only list' % hole
    runs_to_remove_from_dcsonly.remove(hole)

DCSOnly_ll           = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt')
DCSOnlyForNewRuns_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/DCSOnly/json_DCSONLY.txt')
DCSOnlyForNewRuns_ll.removeRuns(runs_to_remove_from_dcsonly)

# Remove runs outside the range [first_run, last_run] since DCS-only
# list includes HI runs, etc.
for ll in (DCSOnly_ll, DCSOnlyForNewRuns_ll):
    ll.removeRuns(xrange(1, first_run))
    ll.removeRuns(xrange(last_run+1, 200000)) # dummy number

Prompt_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-%i_7TeV_PromptReco_Collisions11_JSON%s.txt'          % (last_prompt_run, prompt_version))
PromptMuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-%i_7TeV_PromptReco_Collisions11_JSON_MuonPhys%s.txt' % (last_prompt_run, prompt_version))

May10_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-%i_7TeV_May10ReReco_Collisions11_JSON_v3.txt'          % last_rereco_run)
May10MuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-%i_7TeV_May10ReReco_Collisions11_JSON_MuonPhys_v3.txt' % last_rereco_run)

def combine(may10_ll, prompt_ll, dcsonly_ll=None):
    prompt_ll = copy.deepcopy(prompt_ll)
    prompt_ll.removeRuns(xrange(first_run, last_rereco_run+1))
    ll = may10_ll | prompt_ll
    if dcsonly_ll is not None:
        dcsonly_ll = copy.deepcopy(dcsonly_ll)
        dcsonly_ll.removeRuns(runs_to_remove_from_dcsonly)
        ll = ll | dcsonly_ll
    return ll

Run2011_ll          = combine(May10_ll,          Prompt_ll)
Run2011MuonsOnly_ll = combine(May10MuonsOnly_ll, PromptMuonsOnly_ll)

Run2011PlusDCSOnly_ll          = combine(May10_ll,          Prompt_ll,          DCSOnly_ll)
Run2011PlusDCSOnlyMuonsOnly_ll = combine(May10MuonsOnly_ll, PromptMuonsOnly_ll, DCSOnly_ll)

Run2010_ll          = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON.txt')
Run2010MuonsOnly_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Apr21ReReco_Collisions10_JSON_MuonPhys.txt')
                                   
all_ll_names = ['DCSOnly', 'DCSOnlyForNewRuns', 'Prompt', 'PromptMuonsOnly', 'May10', 'May10MuonsOnly', 'Run2011', 'Run2011MuonsOnly', 'Run2011PlusDCSOnly', 'Run2011PlusDCSOnlyMuonsOnly', 'Run2010', 'Run2010MuonsOnly']

def all_lls():
    return [(x, eval(x + '_ll')) for x in all_ll_names]

for base_name, ll in all_lls():
    exec '%s_ll = ll - to_remove' % base_name
    exec '%s = for_cmssw(%s_ll)' % (base_name, base_name)

if __name__ == '__main__':
    import sys
    if 'write' in sys.argv:
        Run2011MuonsOnly_ll.writeJSON('Run2011MuonsOnly.json')
    elif 'write_all' in sys.argv:
        for base_name, ll in all_lls():
            ll.writeJSON('zp2mu_goodlumis_%s.json' % base_name)
