### Usage:

```sh
bash run_nextflow.sh nextflow.config
```

### Description:
A three step nextflow pipeline that aligns the fastqs to a reference genome and then sorts and generates an index for the aligned bam.

### Setup requirements:
Make sure you configure your AWS CLI installation with appropriate access key id and secret access key. The following example shows sample values. Replace them with your own values.

```sh
$ aws configure
AWS Access Key ID [None]: AKIAIOSFODNN7EXAMPLE
AWS Secret Access Key [None]: *************************RfiCYEXAMPLEKEY
Default region name [None]: us-west-1
Default output format [None]: json
```
