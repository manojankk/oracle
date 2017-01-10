-- DBA Queries

-- To find tables pace
select * from user_segments; or  select * from USER_TABLESPACES;

-- lob
select * from user_lob;

 Find Used and Free Space - User level
SELECT 
USED.TABLESPACE_NAME, 
USED.USED_BYTES AS "USED SPACE(IN GB)",  
FREE.FREE_BYTES AS "FREE SPACE(IN GB)"
FROM
(SELECT TABLESPACE_NAME,TO_CHAR(SUM(NVL(BYTES,0))/1024/1024/1024, '99,999,990.99') AS USED_BYTES FROM USER_SEGMENTS GROUP BY TABLESPACE_NAME) USED
INNER JOIN
(SELECT TABLESPACE_NAME,TO_CHAR(SUM(NVL(BYTES,0))/1024/1024/1024, '99,999,990.99') AS FREE_BYTES FROM  USER_FREE_SPACE GROUP BY TABLESPACE_NAME) FREE
ON (USED.TABLESPACE_NAME = FREE.TABLESPACE_NAME);


-- Find table level memory
select segment_name as table_name,bytes/1024/1024 space_in_mb
 from user_segments
 where segment_type='TABLE'  -------------------------------> object_type
 and segment_name in
 ('FKTA_D_CORRESPONDANT');  -- table_name here....

-- All SCHEMA and TABLESPACE
select * from all_tables;
or
select * from dba_tables;

-- to get ddl
SELECT dbms_metadata.get_ddl('TABLE','WORKFLOW_WORKITEM','LC_CMR_EIM_DEV_01') FROM DUAL;

-- to get db_name, schema, user etc
select 
SYS_CONTEXT('userenv','db_name') as DB_NAME,
SYS_CONTEXT('userenv','session_user') as USERID,
SYS_CONTEXT('userenv','CURRENT_SCHEMA') as schema,
--SYS_CONTEXT('userenv','INSTANCE_NAME') as INSTANCE_NAME,
SYS_CONTEXT('userenv','SERVICE_NAME') as SERVICE_NAME,
SYS_CONTEXT('userenv','IP_ADDRESS') as IP_ADDRESS
from dual;

select * from database_properties; -- 
select * from v$version; -- version 

-- Drop all objects
select 'drop ' || object_type || ' ' || object_name ||
decode(object_type,'TABLE',' cascade constraints','') || ';'  from user_objects
where object_type in 
(
'TABLE',
'SEQUENCE',
'PROCEDURE',
'PACKAGE',
'FUNCTION',
'VIEW'
)
;

-- Find all lob/clob objects,name
SELECT distinct 'purge table ' || table_name || ';' FROM USER_TAB_COLS 
where data_type in ('CLOB','RAW','BLOB');

-- SYS_CONTEXT -- 
ACTION
AUDITED_CURSORID
AUTHENTICATED_IDENTITY
AUTHENTICATION_DATA
AUTHENTICATION_METHOD
BG_JOB_ID
CLIENT_IDENTIFIER
CLIENT_INFO
CURRENT_BIND
CURRENT_SCHEMA
CURRENT_SCHEMAID
CURRENT_SQL
CURRENT_SQLn
CURRENT_SQL_LENGTH
DB_DOMAIN
DB_NAME
DB_UNIQUE_NAME
ENTRYID
ENTERPRISE_IDENTITY
FG_JOB_ID
GLOBAL_CONTEXT_MEMORY
GLOBAL_UID
HOST
IDENTIFICATION_TYPE
INSTANCE
INSTANCE_NAME
IP_ADDRESS
ISDBA
LANG
LANGUAGE
MODULE
NETWORK_PROTOCOL
NLS_CALENDAR
NLS_CURRENCY
NLS_DATE_FORMAT
NLS_DATE_LANGUAGE
NLS_SORT
NLS_TERRITORY
OS_USER
POLICY_INVOKER
PROXY_ENTERPRISE_IDENTITY
PROXY_GLOBAL_UID
PROXY_USER
PROXY_USERID
SERVER_HOST
SERVICE_NAME
SESSION_USER
SESSION_USERID
SESSIONID
SID
STATEMENTID
TERMINAL
