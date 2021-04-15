import argparse
from datetime import datetime
from augur.utils import read_metadata
import json


#class applying SINGLETON Design Pattern
class Recency:
    __instance_recency = None   #instance Variable

    #create instace of class only if previously dont exist else create new
    def __init__(self):
        if Recency.__instance_recency is None:
            Recency.__instance_recency = self
         
        
    @staticmethod
    def get_recency(date_str, ref_date):
        date_submitted = datetime.strptime(date_str, '%Y-%m-%d').toordinal()
        ref_day = ref_date.toordinal()

        delta_days = ref_day - date_submitted
        if delta_days<=0:
            __instance_recency = 'New'
        elif delta_days<3:
            __instance_recency= '1-2 days ago'
        elif delta_days<8:
            __instance_recency= '3-7 days ago'
        elif delta_days<15:
            __instance_recency= 'One week ago'
        elif delta_days<31:
            __instance_recency= 'One month ago'
        elif delta_days>=31:
            __instance_recency= 'Older'

        if not Recency.__instance_recency:
                Recency()         
        return Recency.__instance_recency        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign each sequence a field that specifies when it was added",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help="metadata file")
    parser.add_argument('--output', type=str, required=True, help="output json")
    args = parser.parse_args()

    meta, columns = read_metadata(args.metadata)

    node_data = {'nodes':{}}
    ref_date = datetime.now()
    

    for strain, d in meta.items():
        if 'date_submitted' in d and d['date_submitted']:       
            node_data['nodes'][strain] = {'recency': Recency.get_recency(d['date_submitted'], ref_date)}
    #calling get_recency method 
    
    with open(args.output, 'wt') as fh:
        json.dump(node_data, fh)
