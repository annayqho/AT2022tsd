import numpy as np
import secrets
from penquins import Kowalski


def logon():
    """ Log onto Kowalski """
    username = secrets.un 
    password = secrets.pwd 
    s = Kowalski(
        protocol='https', host='kowalski.caltech.edu', port=443,
            verbose=False, username=username, password=password)
    return s


def get_pos(s, name):
    """ Take median position from alerts """
    det_alerts = get_dets(s, name)
    ras = [det['candidate']['ra'] for det in det_alerts]
    decs = [det['candidate']['dec'] for det in det_alerts]
    ra = np.median(ras)
    dec = np.median(decs)
    return ra,dec


def get_ps1(s, name):
    q = {"query_type": "find_one",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {
                     'objectId': {'$eq': name},
                     'candidate.isdiffpos': {'$in': ['1', 't']},
             },
             "projection": {
                     "_id": 0,
                     "candidate.distpsnr1": 1,
                     "candidate.sgmag1": 1,
             }
         }
         }
    query_result = s.query(query=q)
    try:
        out = query_result['data']['candidate']['distpsnr1']
        out2 = query_result['data']['candidate']['sgmag1']
        return out, out2
    except:
        return []


def get_dets(s, name):
    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {
                     'objectId': {'$eq': name},
                     'candidate.isdiffpos': {'$in': ['1', 't']},
             },
             "projection": {
                     "_id": 0,
                     "candidate.jd": 1,
                     "candidate.magpsf": 1,
                     "candidate.sigmapsf": 1,
                     "candidate.fid": 1,
                     "candidate.programid": 1,
                     "candidate.field": 1,
                     "candidate.ra": 1,
                     "candidate.dec": 1,
             }
         }
         }
    query_result = s.query(query=q)
    try:
        out = query_result['data']
        return out
    except:
        return []


def get_prv_dets(s, name):
    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts_aux",
             "filter": {
                     '_id': {'$eq': name},
             },
             "projection": {
                     "_id": 0,
                     "prv_candidates": 1,
             }
         }
         }
    query_result = s.query(query=q)
    out = query_result['data'][0]['prv_candidates']
    return out


def bin_lc(x,y,ey,bin_size):
    """ Bin a light curve into windows of size bin_size
    Return the binned light curve
    """
    x_binned = []
    y_binned = []
    ey_binned = []

    for ii,x_val in enumerate(x):
        choose = np.abs(x-x_val)<bin_size
        if sum(choose)==1:
            x_binned.append(x_val)
            y_binned.append(y[choose][0])
            ey_binned.append(ey[choose][0])
        elif sum(choose)>1:
            mean,wsum = np.average(
                y[choose], weights=1/ey[choose]**2, returned=True)
            efmean = np.sqrt(1/wsum)
            x_binned.append(np.average(x[choose]))
            y_binned.append(mean)
            ey_binned.append(efmean)

    x_binned = np.array(x_binned)
    y_binned = np.array(y_binned)
    ey_binned = np.array(ey_binned)

    return x_binned,y_binned,ey_binned



def get_lc(s, name):
    """ Retrieve LC for object """
    det_alerts = get_dets(s, name)

    # If the detection resulted in an alert
    isalert = []

    # Detections
    jd = []
    mag = []
    emag = []
    filt = []
    program = []
    fields = []

    # Non-detections
    limjds = []
    limmags = []
    limfilts =[]
    limprogram = []

    for det in det_alerts:
        cand = det['candidate']
        current_jd = cand['jd']
        current_mag = cand['magpsf']
        current_emag = cand['sigmapsf']
        current_filter = cand['fid']
        current_prog = cand['programid']
        current_field = cand['field']

        if current_jd not in jd:
            isalert.append(True)
            jd.append(current_jd)
            mag.append(current_mag)
            emag.append(current_emag)
            filt.append(current_filter)
            program.append(current_prog)
            fields.append(current_field)

    det_prv = get_prv_dets(s, name)
    for prv_cand in det_prv:
        if 'magpsf' in prv_cand.keys():
            if prv_cand['jd'] not in jd:
                if np.logical_or(prv_cand['isdiffpos']=='t', prv_cand['isdiffpos']=='1'):
                    isalert.append(False)
                    jd.append(prv_cand['jd'])
                    mag.append(prv_cand['magpsf'])
                    emag.append(prv_cand['sigmapsf'])
                    filt.append(prv_cand['fid'])
                    program.append(prv_cand['pid'])
                    fields.append(prv_cand['field'])
        else:
            limjds.append(prv_cand['jd'])
            limmags.append(prv_cand['diffmaglim'])
            limfilts.append(prv_cand['fid'])
            limprogram.append(prv_cand['pid'])
    isalert = np.array(isalert)
    jd = np.array(jd)
    mag = np.array(mag)
    emag = np.array(emag)
    filt = np.array(filt)
    program = np.array(program)
    fields = np.array(fields)
    limjds = np.array(limjds)
    limmags = np.array(limmags)
    limfilts = np.array(limfilts)
    limprogram = np.array(limprogram)

    # Sort in order of jd
    order = np.argsort(jd)
    isalert = isalert[order]
    jd = jd[order]
    mag = mag[order]
    emag = emag[order]
    filt = filt[order]
    fields = fields[order]
    program = program[order]

    order = np.argsort(limjds)
    limjds = limjds[order]
    limmags = limmags[order]
    limfilts = limfilts[order]
    limprogram = limprogram[order]

    return isalert,jd,mag,emag,filt,program,limjds,limmags,limfilts,limprogram
