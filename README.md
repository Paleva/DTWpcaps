# Overview
This project implements network analysis of traffic between a target and an attacker by using pcaps from respective machines. The analysis consists of using three algorithms to detect similarities between the captures. In this particular case the exploit used is [CVE-2020-10564](https://nvd.nist.gov/vuln/detail/CVE-2020-10564) with [POC](https://github.com/beerpwn/CVE/blob/master/WP-File-Upload_disclosure_report/CVE-2020-10564_exploit.py).


# Usage

Create python venv and install the required packages in requirements.txt

``` bash
python3 analysis.py --attacker <attacker_pcap> --target <target_pcap>
```

# Result 
This produces multiple graphs that show what each algorithm produces in regards to the similarity and can help identify anomalies in traffic if there's noise or auxilary traffic.